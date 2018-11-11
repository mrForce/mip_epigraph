from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACProtein
from Bio.Alphabet import IUPAC


def clean_proteins(protein_location):
    with open(protein_location, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            record.seq = Seq(''.join(list(filter(lambda x: x in IUPACProtein.letters, str(record.seq)))), IUPAC.protein)
            print(record.format('fasta'))


clean_proteins('HIV_small_joined.fasta')
