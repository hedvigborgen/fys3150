from formats import Table

tab = [[]]
with open('table.txt', 'r') as infile:
    tab[0] = [word.replace('<space>', ' ') for word in infile.readline().split()]
    for line in infile:
        tab.append([word.replace('<space>', ' ') for word in line.split()])
tab = Table(tab)
tab.write()
tab.latex(filename='table.tex', complete=True)