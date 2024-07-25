import gzip
import sys

def create_fake_fastq(sample_name):
    fake_fastq_content = f"""@{sample_name}_SEQ_ID
N
+
#
"""
    file_name = f'{sample_name}_fq_end2.fastq.gz'
    with gzip.open(file_name, 'wt') as f:
        f.write(fake_fastq_content)
    print(f"Created {file_name}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script_name.py <sample_name1> <sample_name2> ...")
    else:
        sample_names = sys.argv[1:]
        for sample_name in sample_names:
            create_fake_fastq(sample_name)
