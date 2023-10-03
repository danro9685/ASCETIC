ASCETIC
=======

| Branch | Status |
| --- | --- |
| master | [![R-CMD-check-bioc](https://github.com/danro9685/ASCETIC/actions/workflows/check-bioc.yml/badge.svg?branch=master)](https://github.com/danro9685/ASCETIC/actions/workflows/check-bioc.yml) |
| development | [![R-CMD-check-bioc](https://github.com/danro9685/ASCETIC/actions/workflows/check-bioc.yml/badge.svg?branch=development)](https://github.com/danro9685/ASCETIC/actions/workflows/check-bioc.yml) |

Cancer development is a stochastic process involving large populations of cells. Random genetic and epigenetic alterations commonly occurring in any cell can occasionally be beneficial to neoplastic ones, thus defining clones characterized by a functional selective advantage. During clonal evolution, certain clones can be positively selected for increased proliferation and survival ability, outgrowing competing cells and this can eventually lead to invasion and metastasis. Throughout such a multi-step stochastic process, cancer cells can acquire over time a set of biological capabilities that sometimes are referred to as hallmarks. Not all variants are involved in their acquisition, but only a relative small subset of them – i.e., the drivers –, while most mutations present in the cancer clones do not increase their fitness – i.e., the passengers.

In response to the constantly increasing availability of cancer omics data and the rapid advancement of data science and machine learning techniques, we introduce a novel framework named ASCETIC (Agony-baSed Cancer EvoluTion InferenCe, https://www.nature.com/articles/s41467-023-41670-3). ASCETIC can extract cancer's evolutionary signatures, which represent recurring paths of driver mutation acquisition associated with distinct disease outcomes. ASCETIC can process sequencing data derived from different technologies, including bulk and single-cell sequencing. The approach goes beyond the traditional focus on individual genetic alterations and explores the broader landscape of genomic alterations and their intricate interactions, with the aim to enhance predictive accuracy and our understanding of cancer evolution's impact on prognosis.

ASCETIC's workflow involves the reconstruction of robust tumor evolution models for individual patients, followed by the integration of these models into a comprehensive cancer-specific evolution model. By leveraging prognostic data through regularized Cox regression, ASCETIC identifies significant evolutionary patterns or signatures associated with patient outcomes.

In summary, ASCETIC represents a powerful tool for uncovering consistent evolutionary patterns in cancer, offering the potential to contribute to a curated catalogue of evolutionary signatures akin to the widely used COSMIC Mutational Signatures database.

Please feel free to contact us if you have any questions regarding our tool at daniele.ramazzotti@unimib.it.
