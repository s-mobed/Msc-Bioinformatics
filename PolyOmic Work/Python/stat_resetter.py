
with open('WT_dyno_stats.csv') as stats:
    for line in stats:
        comma = line.replace(' ', ',')
        comma_list = comma.split(',')
        new_line = []

        for entry in comma_list:
            if len(entry) > 0 and entry != '\n':
                new_line.append(entry)


        final_line = '\t'.join(new_line) + '\n'
        print(final_line)


