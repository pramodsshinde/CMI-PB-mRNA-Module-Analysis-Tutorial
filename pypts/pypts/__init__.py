import pandas as pd 

def load_api_data(url):
    """
    Loading data using the API.
    """

    import urllib.request, json, ssl 

    gcontext = ssl.SSLContext()  
    with urllib.request.urlopen(url, context=gcontext) as url:
        data = url.read()
        data = pd.read_json(data)
        return(data)

def load_subject_sample_map():
    """
    Loading the subject to samples merged table.
    """
    subjects = load_api_data('https://staging.cmi-pb.org:443/db/subject')
    samples = load_api_data('https://staging.cmi-pb.org:443/db/sample')
    master = subjects.merge(samples, on='subject_id')
    return(master)


def load_innatedb_curated_genes():
    
    """
    Load data which was downloaded from the InnateDB website. 
    """
    
    immune_genes = pd.read_table('../data/innatedb_curated_genes.txt', index_col=0)

    # By digging around the InnateDB I found out that the species ID for humans is 9606
    immune_genes = immune_genes[immune_genes['Species'] == 9606]
    immune_genes.drop(['Annotation', 'PubMED ID'], inplace=True, axis=1)
    immune_genes.drop_duplicates(subset=['Gene Symbol'], inplace=True)
    immune_genes = immune_genes['Gene Symbol'].squeeze().values
    
    return(immune_genes)
    
def preboost_conversion(x):
    if x > 0:
        return(x)
    else:
        return('preboost') 
