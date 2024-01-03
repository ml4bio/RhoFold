from Bio.Blast import NCBIXML

# 解析BLASTN结果文件
with open("blast_result.xml", "r") as result_file:
    blast_record = NCBIXML.read(result_file)

# 提取第一条比对记录的相似度百分比
similarity = blast_record.alignments[0].hsps[0].identities / blast_record.alignments[0].hsps[0].align_length * 100
print("Similarity: {:.2f}%".format(similarity))