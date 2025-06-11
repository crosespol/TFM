import json
import pandas as pd
from pathlib import Path

# Paths
csv_path = Path("/home/carlesroses/PHROGs/accessions.csv")
template_path = Path("/home/carlesroses/PHROGs/CGA-MD/run/config/config.json")
output_dir = Path("/home/carlesroses/PHROGs/CGA-MD/run/config/")

# Carrega plantilla
with open(template_path) as f:
    template = json.load(f)

# Llegeix CSV (separat per comes)
df = pd.read_csv(csv_path, dtype=str)

# Genera un .json per cada fila
for _, row in df.iterrows():
    genome_id = row["genome_id"]
    cluster_size = row["cluster_size"]
    protein_id = row["protein_id"]

    config = template.copy()
    config["input_records"] = [protein_id]
    config["output_parameters"]["folder_name"] = f"output_cluster{cluster_size}_{genome_id}_{protein_id}"

    json_path = output_dir / f"config_{protein_id}.json"
    with open(json_path, "w") as f:
        json.dump(config, f, indent=4)

    print(f"[✓] Creat: {json_path}")
