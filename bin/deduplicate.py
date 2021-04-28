from Bio import SeqIO

seen = []
out_handle = open("/Users/rmcolq/Work/git/grapevine-global-tree/test/test2.fa", "w")
record_dict = SeqIO.parse("/Users/rmcolq/Work/git/grapevine-global-tree/test/test.fa", "fasta")
for record in record_dict:
    if record.id in seen:
        continue
    SeqIO.write(record, out_handle, "fasta-2line")
    seen.append(record.id)
out_handle.close()
