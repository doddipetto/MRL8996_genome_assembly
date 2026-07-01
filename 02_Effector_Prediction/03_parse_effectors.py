import os

# File di input e output
deeptmhmm_file = "biolib_results/predicted_topologies.3line"
fasta_in = "SignalP_EffectorP_intersect.fasta"
fasta_out = "FoxMRL8996_Final_Effectorome.fasta"
ids_out = "FoxMRL8996_Final_Effectorome_IDs.txt"

print("Analisi topologie DeepTMHMM in corso...")

valid_ids = set()
with open(deeptmhmm_file, 'r') as f:
    for line in f:
        if line.startswith(">"):
            # Esempio: >g1.t1 | type=SP+TM
            header = line.strip().split(" | ")
            prot_id = header[0][1:] 
            topology = header[1] if len(header) > 1 else ""
            
            # Scartiamo chiunque abbia 'TM' nella topologia
            if "TM" not in topology:
                valid_ids.add(prot_id)

print(f"Sopravvissuti al filtro di membrana: {len(valid_ids)} candidati.")

# Generazione del FASTA finale
print("Estrazione sequenze definitive...")
kept_count = 0
with open(fasta_in, 'r') as f_in, open(fasta_out, 'w') as f_out, open(ids_out, 'w') as id_out:
    keep = False
    for line in f_in:
        if line.startswith(">"):
            seq_id = line.strip().split()[0][1:] 
            if seq_id in valid_ids:
                keep = True
                kept_count += 1
                f_out.write(line)
                id_out.write(seq_id + "\n")
            else:
                keep = False
        elif keep:
            f_out.write(line)

print("\n--- EFFETTOROMA FoxMRL8996 COMPLETATO ---")
print(f"Proteine totali: {kept_count}")
print(f"FASTA salvato in: {fasta_out}")
print(f"Lista ID salvata in: {ids_out}")
