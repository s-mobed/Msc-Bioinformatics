import argparse
import sys

interface = argparse.ArgumentParser(prog='VCF file parser',
                                    description='This program takes a VCF, GFF, and fasta file to generate'
                                                'a polished final table of the SNPs stored in the VCF file'
                                                'and the effect of those mutations. A bar plot for the counts'
                                                'of non-coding, non-synonymous, and synonymous mutations.')
interface.add_argument('filename', nargs='*', help='Input files can only be either GFF/VCF/fasta formats')

interface_args = interface.parse_args()
file_list_str = ''
file_list_ext = [element.split('.')[1] for element in interface_args.filename]
print(file_list_ext)
if len(interface_args.filename) == 0:
    interface.print_help()
    print('='*42)
    sys.exit('Solution: Please consider inputting files' + '\n' + '='*42)

elif len(file_list_ext) != len(set(file_list_ext)):
    print('=' * 42)
    sys.exit('Duplicate file formats detected, only insert one of each file types' + '\n' + '=' * 42)

elif set(file_list_ext) != {'gff', 'fasta', 'vcf'}:
    print('=' * 42)
    missing = [ext for ext in ['fasta', 'gff', 'vcf'] if ext not in file_list_ext]
    for missing in missing:
        print(f'Missing {missing} file')
    sys.exit(f'Exit due to missing above files' + '\n' + '=' * 42)

else:
    for file in interface_args.filename:
        if file.endswith(('.vcf', 'vcf.gz')):
            vcf_data = file
        elif file.endswith('.gff'):
            gff_data = file
        elif file.endswith('.fasta'):
            fasta_data = file
        else:
            interface.print_help()
            print('=' * 42)
            sys.exit(f'{file} is not an acceptable file format' + '\n' + '=' * 42)

        file_list_str += file.split("/")[-1] + ', '

