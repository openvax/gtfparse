from data import data_path
from gtfparse import read_gtf


ENSEMBL_GTF_PATH = data_path("ensembl_grch37.head.gtf")

EXPECTED_FEATURES = set([
    "gene",
    "transcript",
    "exon",
    "CDS",
    "UTR",
    "start_codon",
    "stop_codon",
])


def test_ensembl_gtf_columns():
    df = read_gtf(ENSEMBL_GTF_PATH)
    features = set(df["feature"])
    assert features == EXPECTED_FEATURES

# first 1000 lines of GTF only contained these genes
EXPECTED_GENE_NAMES = {
    'FAM41C', 'CICP27', 'RNU6-1100P', 'NOC2L', 'AP006222.1',
    'LINC01128', 'RP4-669L17.1', 'RP11-206L10.2', 'PLEKHN1',
    'WBP1LP7', 'RP5-857K21.1', 'RP5-857K21.5', 'RNU6-1199P',
    'RP11-206L10.10', 'RP11-54O7.16', 'CICP7', 'AL627309.1',
    'RP5-857K21.11', 'DDX11L1', 'RP5-857K21.3', 'RP11-34P13.7',
    'AL669831.1', 'MTATP6P1', 'CICP3', 'WBP1LP6', 'LINC00115',
    'hsa-mir-6723', 'RP5-857K21.7', 'SAMD11', 'RP11-206L10.5',
    'RP11-34P13.8', 'RP11-206L10.9', 'RP11-34P13.15', 'TUBB8P11',
    'MTATP8P1', 'RP4-669L17.8', 'RP11-206L10.1', 'RP11-34P13.13',
    'RP11-206L10.3', 'RP11-206L10.4', 'RP11-54O7.3', 'RP5-857K21.2',
    'OR4F5', 'MTND1P23', 'AL645608.1', 'RP11-34P13.16', 'RP11-34P13.14',
    'AP006222.2', 'OR4F29', 'RP4-669L17.4', 'AL732372.1', 'OR4G4P',
    'MTND2P28', 'OR4F16', 'KLHL17', 'FAM138A', 'OR4G11P', 'FAM87B',
    'RP5-857K21.15', 'AL645608.2', 'RP11-206L10.8', 'RP5-857K21.4',
    'MIR1302-10', 'RP11-54O7.2', 'RP4-669L17.10', 'RP11-54O7.1',
    'RP11-34P13.9', 'WASH7P', 'RP4-669L17.2'
}

def test_ensembl_gtf_gene_names():
    df = read_gtf(ENSEMBL_GTF_PATH)
    gene_names = set(df["gene_name"])
    assert gene_names == EXPECTED_GENE_NAMES, \
        "Wrong gene names: %s, missing %s and unexpected %s" % (
            gene_names,
            EXPECTED_GENE_NAMES.difference(gene_names),
            gene_names.difference(EXPECTED_GENE_NAMES)
        )

def test_ensembl_gtf_gene_names_with_usecols():
    df = read_gtf(ENSEMBL_GTF_PATH, usecols=["gene_name"])
    gene_names = set(df["gene_name"])
    assert gene_names == EXPECTED_GENE_NAMES, \
        "Wrong gene names: %s, missing %s and unexpected %s" % (
            gene_names,
            EXPECTED_GENE_NAMES.difference(gene_names),
            gene_names.difference(EXPECTED_GENE_NAMES)
        )

def test_ensembl_gtf_gene_names_zip():
    df = read_gtf(ENSEMBL_GTF_PATH + ".gz")
    gene_names = set(df["gene_name"])
    assert gene_names == EXPECTED_GENE_NAMES, \
        "Wrong gene names: %s, missing %s and unexpected %s" % (
            gene_names,
            EXPECTED_GENE_NAMES.difference(gene_names),
            gene_names.difference(EXPECTED_GENE_NAMES)
        )

def test_ensembl_gtf_gene_names_with_usecols_gzip():
    df = read_gtf(ENSEMBL_GTF_PATH + ".gz", usecols=["gene_name"])
    gene_names = set(df["gene_name"])
    assert gene_names == EXPECTED_GENE_NAMES, \
        "Wrong gene names: %s, missing %s and unexpected %s" % (
            gene_names,
            EXPECTED_GENE_NAMES.difference(gene_names),
            gene_names.difference(EXPECTED_GENE_NAMES)
        )
