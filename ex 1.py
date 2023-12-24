def square_matrix_elements(input_file, output_file):
    with open(input_file, 'r') as file:
        matrix_lines = file.readlines()
    matrix = []
    for line in matrix_lines:
        row = [int(x) for x in line.strip().split()]
        squared_row = [x ** 2 for x in row]
        matrix.append(squared_row)

    with open(output_file, 'w') as file:
        for row in matrix:
            line = ' '.join(str(x) for x in row)
            file.write(line + '\n')

square_matrix_elements("input.txt", "out.dat")
