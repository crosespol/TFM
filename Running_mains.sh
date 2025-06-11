#!/bin/bash

ROOT_DIR="/home/carlesroses/PHROGs/CGA-MD/run"
CONFIG_DIR="$ROOT_DIR/config"
LOG_DIR="$ROOT_DIR/logs"
SCRIPT_PATH="/home/carlesroses/PHROGs/CGA-MD/scripts/main.py"

mkdir -p "$LOG_DIR"

for json_file in "$CONFIG_DIR"/config_*.json; do
    name=$(basename "$json_file" .json)
    log_file="$LOG_DIR/${name}.log"

    echo "[RUNNING] $name"

    # Copiem aquest fitxer com config.json
    cp "$json_file" "$CONFIG_DIR/config.json"

    # Executem amb el root directory
    python3 "$SCRIPT_PATH" "$ROOT_DIR" > "$log_file" 2>&1

    echo "[✓] Finalitzat: $name"
done

echo "[✓] Totes les execucions s'han completat."
