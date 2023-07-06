
import collections

mod_count = 0
high_count = 0
ids_list = []
with open('sifted_AmpB_compliment.vcf') as file:
    with open('sifted_id.txt', 'w') as ids:
        with open('high_sifted_id.txt', 'w') as high_output:
            with open('het_sifted_id.txt', 'w') as het_output:
                for line in file:
                    line_list = line.split(sep='\t')
                    if line_list[0] != 'LmxM.00':
                        # transcript ids
                        ids.write(line_list[6] + '\n')
                        ids_list.append(line_list[6])

                        # High impacts snps
                        if line_list[4] == 'HIGH':
                            high_count += 1
                            print(line)
                            high_output.write(line_list[6] + '\n')

                        # Moderate impact snps
                        if line_list[4] == 'MODERATE':
                            mod_count += 1

                        # Het snps
                        elif line_list[9].startswith('0/1'):
                            print(line)
                            het_output.write(line_list[6] + '\n')

counts = collections.Counter(ids_list)

# Sort the counts in descending order and print them
sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)
for element, count in sorted_counts:
    if count >1 :
        print(element)