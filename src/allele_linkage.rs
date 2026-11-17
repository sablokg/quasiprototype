use petgraph::graph::{NodeIndex, UnGraph};
use petgraph::unionfind::UnionFind;
use rust_htslib::bam::{self, Read as BamRead};
use rust_htslib::bcf::{self, Read as BcfRead};
use rust_htslib::faidx;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Write};

/*
Gaurav Sablok
gsablok@proton.me
*/

#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct AlleleNode {
    pub chrom: String,
    pub pos: u64, // 1-based genomic position
    pub base: char,
}

#[derive(Debug, Clone, Copy)]
pub struct EdgeWeight {
    pub phi: f64,
    pub n11: usize,
}

pub type ReadObs = HashMap<Vec<u8>, BTreeMap<(String, u64), char>>;

/// 1. Stream BAM + VCF to extract read-to-allele observations
pub fn load_reads_from_bam(
    bam_path: &str,
    vcf_path: &str,
    min_mapq: u8,
    min_baseq: u8,
) -> Result<ReadObs, Box<dyn Error>> {
    let mut bam = bam::IndexedReader::from_path(bam_path)?;
    let mut vcf = bcf::Reader::from_path(vcf_path)?;
    let header = vcf.header().clone();

    let mut reads: ReadObs = HashMap::new();

    for record_res in vcf.records() {
        let record = record_res?;
        let alleles = record.alleles();

        // SNPs only
        if alleles.len() < 2 || alleles.iter().any(|a| a.len() != 1) {
            continue;
        }
        let valid_bases: HashSet<u8> = alleles.iter().map(|a| a[0].to_ascii_uppercase()).collect();

        let rid = match record.rid() {
            Some(r) => r,
            None => continue,
        };
        let chrom = String::from_utf8_lossy(header.rid2name(rid)?).into_owned();
        let pos0: i64 = record.pos();
        let pos1 = (pos0 + 1) as u64;

        // Contig missing from BAM -> skip this site
        if bam.fetch((chrom.as_str(), pos0, pos0 + 1)).is_err() {
            continue;
        }

        // Per-site observations; None marks mates that disagree
        let mut site_obs: HashMap<Vec<u8>, Option<char>> = HashMap::new();

        for rec_res in bam.records() {
            let rec = rec_res?;
            if rec.is_unmapped()
                || rec.is_secondary()
                || rec.is_supplementary()
                || rec.is_duplicate()
                || rec.is_quality_check_failed()
                || rec.mapq() < min_mapq
            {
                continue;
            }

            // Query position aligned to the SNP (None for deletion / skip / not covered)
            let qpos = match rec.cigar().read_pos(pos0 as u32, false, false)? {
                Some(q) => q as usize,
                None => continue,
            };

            let qual = rec.qual()[qpos];
            let base = rec.seq().as_bytes()[qpos].to_ascii_uppercase();
            if qual < min_baseq || !valid_bases.contains(&base) {
                continue;
            }

            let base = base as char;
            site_obs
                .entry(rec.qname().to_vec())
                .and_modify(|e| {
                    if *e != Some(base) {
                        *e = None; // overlapping mates conflict
                    }
                })
                .or_insert(Some(base));
        }

        for (qname, base) in site_obs {
            if let Some(b) = base {
                reads
                    .entry(qname)
                    .or_default()
                    .insert((chrom.clone(), pos1), b);
            }
        }
    }

    Ok(reads)
}

/// 2. Build the allele-linkage graph
pub fn build_allele_graph(
    reads: &ReadObs,
    min_phi: f64,
    min_n11: usize,
) -> UnGraph<AlleleNode, EdgeWeight> {
    type Site = (String, u64);
    let mut pair_counts: BTreeMap<(Site, Site), BTreeMap<(char, char), usize>> = BTreeMap::new();

    for obs in reads.values() {
        let sites: Vec<(&Site, &char)> = obs.iter().collect();
        for i in 0..sites.len() {
            for j in (i + 1)..sites.len() {
                let (s1, &b1) = sites[i];
                let (s2, &b2) = sites[j];
                *pair_counts
                    .entry((s1.clone(), s2.clone()))
                    .or_default()
                    .entry((b1, b2))
                    .or_insert(0) += 1;
            }
        }
    }

    let mut graph = UnGraph::<AlleleNode, EdgeWeight>::new_undirected();
    let mut node_indices: HashMap<AlleleNode, NodeIndex> = HashMap::new();

    for ((s1, s2), counts) in &pair_counts {
        let n: usize = counts.values().sum();
        if n < min_n11 {
            continue;
        }

        // Marginal counts for each allele at each site
        let mut m1: HashMap<char, usize> = HashMap::new();
        let mut m2: HashMap<char, usize> = HashMap::new();
        for (&(a, b), &c) in counts {
            *m1.entry(a).or_insert(0) += c;
            *m2.entry(b).or_insert(0) += c;
        }

        let nf = n as f64;
        for (&(u, v), &n11) in counts {
            if n11 < min_n11 {
                continue;
            }
            let n1 = m1[&u] as f64;
            let n2 = m2[&v] as f64;

            let denom = (n1 * n2 * (nf - n1) * (nf - n2)).sqrt();
            if denom == 0.0 {
                continue; // monomorphic in this read set
            }

            let phi = (n11 as f64 * nf - n1 * n2) / denom;
            if phi < min_phi {
                continue;
            }

            let node_u = AlleleNode {
                chrom: s1.0.clone(),
                pos: s1.1,
                base: u,
            };
            let node_v = AlleleNode {
                chrom: s2.0.clone(),
                pos: s2.1,
                base: v,
            };

            let idx_u = *node_indices
                .entry(node_u.clone())
                .or_insert_with(|| graph.add_node(node_u));
            let idx_v = *node_indices
                .entry(node_v.clone())
                .or_insert_with(|| graph.add_node(node_v));

            graph.add_edge(idx_u, idx_v, EdgeWeight { phi, n11 });
        }
    }

    graph
}

/// 3. Export connected components as haplotype sequences
pub fn export_haplotypes_to_fasta(
    graph: &UnGraph<AlleleNode, EdgeWeight>,
    ref_fa_path: &str,
    out_fa_path: &str,
) -> Result<(), Box<dyn Error>> {
    let fa_reader = faidx::Reader::from_path(ref_fa_path)?;
    let mut writer = BufWriter::new(File::create(out_fa_path)?);

    // Connected components via union-find
    let mut uf = UnionFind::<usize>::new(graph.node_count());
    for e in graph.raw_edges() {
        uf.union(e.source().index(), e.target().index());
    }

    let mut comps: BTreeMap<usize, Vec<NodeIndex>> = BTreeMap::new();
    for n in graph.node_indices() {
        comps.entry(uf.find(n.index())).or_default().push(n);
    }

    // Edge support per component
    let mut support: HashMap<usize, (usize, usize)> = HashMap::new(); // root -> (sum n11, edges)
    for e in graph.raw_edges() {
        let root = uf.find(e.source().index());
        let s = support.entry(root).or_insert((0, 0));
        s.0 += e.weight.n11;
        s.1 += 1;
    }

    let mut ordered: Vec<(usize, Vec<NodeIndex>)> = comps.into_iter().collect();
    ordered.sort_by_key(|(root, c)| (std::cmp::Reverse(c.len()), *root));

    let mut ref_cache: HashMap<String, Vec<u8>> = HashMap::new();

    for (idx, (root, comp)) in ordered.iter().enumerate() {
        let (total, edges) = support.get(root).copied().unwrap_or((0, 0));
        let mean_support = total as f64 / edges.max(1) as f64;

        // chrom -> pos -> set of alleles seen in this component
        let mut by_chrom: BTreeMap<&str, BTreeMap<u64, Vec<char>>> = BTreeMap::new();
        for &ni in comp {
            let node = &graph[ni];
            by_chrom
                .entry(node.chrom.as_str())
                .or_default()
                .entry(node.pos)
                .or_default()
                .push(node.base);
        }

        for (chrom, sites) in by_chrom {
            if !ref_cache.contains_key(chrom) {
                let len = fa_reader.fetch_seq_len(chrom);
                if len == 0 {
                    continue;
                }
                let seq = fa_reader.fetch_seq(chrom, 0, len as usize - 1)?;
                ref_cache.insert(
                    chrom.to_string(),
                    seq.iter().map(|b| b.to_ascii_uppercase()).collect(),
                );
            }
            let mut seq = ref_cache[chrom].clone();

            let mut applied = 0usize;
            let mut conflicts = 0usize;
            for (pos, bases) in &sites {
                let mut distinct = bases.clone();
                distinct.sort_unstable();
                distinct.dedup();
                if distinct.len() > 1 {
                    conflicts += 1; // both alleles linked: ambiguous, keep reference
                    continue;
                }
                let zero_idx = (*pos as usize) - 1;
                if zero_idx < seq.len() {
                    seq[zero_idx] = distinct[0] as u8;
                    applied += 1;
                }
            }

            writeln!(
                writer,
                ">Haplotype_{}_{}|sites={}|conflicts={}|mean_support={:.1}",
                idx + 1,
                chrom,
                applied,
                conflicts,
                mean_support
            )?;
            for line in seq.chunks(80) {
                writer.write_all(line)?;
                writer.write_all(b"\n")?;
            }
        }
    }

    writer.flush()?;
    Ok(())
}
