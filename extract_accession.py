from Bio import Entrez, SeqIO

Entrez.email = "EMAIL@ADRESS"
Entrez.api_key = "API_KEY"
import time

def get_protein_id(genome_id, cds_end):
    try:
        handle = Entrez.efetch(db="nuccore", id=genome_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        

        for feature in record.features:
            if feature.type == "CDS":
                end = int(feature.location.end)
                if end == int(cds_end):
                    return feature.qualifiers.get("protein_id", ["NA"])[0]
        return "NOT_FOUND"
    except Exception as e:
        return f"ERROR: {e}"

# Example usage
test_list = [
("LC727894" , 13737),
("MN270273" , 8325),
("NC_048045" , 80605),
("MH431937" , 21537),
("MK511051" , 43165),
("KM389222" , 811),
("MK448620" , 4105),
("MK448873" , 3213),
("MK448946" , 2753),
("KT337359" , 29402),
("MN534318" , 2866),
("MH431931" , 27384),
("MK448864" , 3738),
("OQ594956" , 12686),
("KY065488" , 3102),
("LN997803" , 24923),
("OK490431" , 3817),
("LC645430" , 11094),
("CP124941" , 26575),
("LR797923" , 15218),
("ON391949" , 3981),
("MK510977" , 26648),
("OP168678" , 25006),
("LC616079" , 6422),
("MN534319" , 2861),
("NC_007805" , 26858),
("LC567822" , 7076),
("AF547987" , 27419),
("MK448767" , 2552),
("NC_005857" , 39145),
("MK448507" , 4511)
]

with open("accessions.csv", "w") as out:
    out.write("genome_id\tend\tprotein_id\n")
    for genome_id, end in test_list:
        prot_id = get_protein_id(genome_id, end)
        print(f"{genome_id}\t{end}\t{prot_id}")
        out.write(f"{genome_id}\t{end}\t{prot_id}\n")
        time.sleep(0.4)
