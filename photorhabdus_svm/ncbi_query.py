import csv
query = ''
number = 0

with open('/Users/desho/Desktop/OG0000017.csv') as file:
    f = csv.reader(file)
    print(f)
    for row in f:
        query = query + row[0] + ' OR '
        number += 1

    print(query)
    print(number)