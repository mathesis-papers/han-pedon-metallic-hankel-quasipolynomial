import argparse
import hashlib
import json
import math
from pathlib import Path


def add(left, right, scale=1):
    result = [0] * max(len(left), len(right))
    for index, value in enumerate(left):
        result[index] += value
    for index, value in enumerate(right):
        result[index] += scale * value
    while len(result) > 1 and result[-1] == 0:
        result.pop()
    return result


def multiply(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for left_index, left_value in enumerate(left):
        if left_value:
            for right_index, right_value in enumerate(right):
                result[left_index + right_index] += left_value * right_value
    return result


def coefficient(polynomial, index):
    return polynomial[index] if 0 <= index < len(polynomial) else 0


def fraction_data(parameter):
    if parameter == 1:
        return [(0, -1, [1, 1]), (0, 1, [1, 1]), (1, -1, [1, 1, -1])]
    linear = add([1, -1], [0] * parameter + [1, -1])
    total = add(linear, [0] + [1] * parameter)
    end = add(total, [0] * (parameter + 1) + [1])
    data = [(parameter - 1, -1, end), (0, 1, [1])]
    for offset in range(parameter - 2):
        data.extend([(parameter - offset - 2, -1, add([1] * (parameter - offset), [0, -1])),
                     (offset, -1, [1] * (offset + 2)), (0, -1, [1, -1])])
    data.append((0, -1, [1]))
    for offset in range(parameter - 1):
        data.extend([(parameter - offset - 2, -1, [1] * (parameter - offset)),
                     (offset, -1, add([1] * (offset + 2), [0, -1])),
                     (0, -1, [1] if offset == parameter - 2 else [1, -1])])
    data.extend([(parameter - 1, 1, end), (parameter, -1, total)])
    return data


def period(parameter):
    return (2 if parameter % 2 == 0 else 4) * parameter * (parameter + 1)


def continuants(parameter, block):
    data = fraction_data(parameter)
    cycles = block * (1 if parameter % 2 == 0 else 2)
    previous, current = [0], [1]
    size, previous_size = 0, 0
    last_order = data[-1][0]
    for position in range(cycles * len(data)):
        order, sign, denominator = data[position % len(data)]
        following = add(multiply(denominator, current), [0] * (last_order + order + 2) + previous, -sign)
        previous, current = current, following
        previous_size, size = size, size + order + 1
        last_order = order
    assert size == block * period(parameter)
    assert len(current) <= size + 1 and len(previous) <= previous_size + 1
    first = [coefficient(current, size - index) for index in range(size + 1)]
    second = [coefficient(previous, previous_size - index) for index in range(previous_size + 1)]
    assert first[-1] == second[-1] == 1
    return first, second, current, [0] * (parameter + 1) + previous


def moments(parameter, count):
    linear = add(add([1, -1], [0] * parameter + [1, -1]), [0] + [1] * parameter, -1)
    values = [1]
    for degree in range(1, count):
        convolution = sum(values[index] * values[degree - 1 - index] for index in range(degree))
        remainder = sum(linear[index] * values[degree - index] for index in range(1, min(degree, len(linear) - 1) + 1))
        values.append(-convolution - remainder)
    return values


def determinant(matrix):
    size = len(matrix)
    if size == 0:
        return 1
    work = [row[:] for row in matrix]
    previous, sign = 1, 1
    for pivot_index in range(size - 1):
        if work[pivot_index][pivot_index] == 0:
            replacement = next((row for row in range(pivot_index + 1, size) if work[row][pivot_index]), None)
            if replacement is None:
                return 0
            work[pivot_index], work[replacement] = work[replacement], work[pivot_index]
            sign = -sign
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = work[row][column] * pivot - work[row][pivot_index] * work[pivot_index][column]
                quotient, remainder = divmod(numerator, previous)
                assert remainder == 0
                work[row][column] = quotient
            work[row][pivot_index] = 0
        previous = pivot
    return sign * work[-1][-1]


def hankel(values, shift, size):
    return determinant([[values[shift + row + column] for column in range(size)] for row in range(size)])


def boundary(first, second, shift_offset, size, dimension):
    normal_size = size + dimension
    reverse_first = [coefficient(first, normal_size - depth) for depth in range(normal_size + 1)]
    reverse_second = [coefficient(second, normal_size - depth) for depth in range(normal_size + 1)]
    matrix = []
    for row in range(dimension):
        output = []
        for column in range(dimension):
            if row < shift_offset:
                value = 0
                for index in range(row + 1):
                    depth = dimension + index - column - row - 1
                    if depth >= 0:
                        value -= coefficient(reverse_first, depth) * coefficient(second, index)
                        value += coefficient(reverse_second, depth) * coefficient(first, index)
            else:
                value = 0
                for index in range(dimension - column):
                    depth = 2 * dimension - 1 - row - column - index
                    value += coefficient(reverse_first, depth) * coefficient(reverse_second, index)
                    value -= coefficient(reverse_second, depth) * coefficient(reverse_first, index)
            output.append(value)
        matrix.append(output)
    return (-1) ** (shift_offset * size) * determinant(matrix)


def differences(values):
    return [right - left for left, right in zip(values, values[1:])]


def golden_polynomial_identities():
    constants = [1, 2, 0, -2, -3, -4, 0, 4]
    slopes = [4, 4, 0, -4, -4, -4, 0, 4]
    endpoints = {}
    identities = []
    for residue in range(8):
        initial_block = (residue + 2 + 7) // 8
        dimension = initial_block * 8 - residue
        degree_bound = dimension * (dimension - 2) + 2
        values = []
        for argument in range(degree_bound + 1):
            block = initial_block + argument
            if block not in endpoints:
                endpoints[block] = continuants(1, block)[:2]
            first, second = endpoints[block]
            value = boundary(first, second, 2, residue + 8 * argument, dimension)
            assert value == constants[residue] + slopes[residue] * argument
            values.append(value)
        identities.append({'residue': residue, 'complementary_dimension': dimension,
                           'degree_bound': degree_bound,
                           'constant': constants[residue], 'slope': slopes[residue],
                           'values_m_0_through_degree_bound': values})
    assert sum(len(row['values_m_0_through_degree_bound']) for row in identities) == 220
    return identities


def modular_newton_identities():
    count = 0
    for degree in range(9):
        denominator = math.lcm(*range(1, max(1, degree) + 1))
        for modulus in range(1, 21):
            step = modulus * denominator
            for argument in range(-4, 5):
                for order in range(degree + 1):
                    first = (math.comb(argument, order) if argument >= 0
                             else (-1) ** order * math.comb(-argument + order - 1, order))
                    shifted = argument + step
                    second = (math.comb(shifted, order) if shifted >= 0
                              else (-1) ** order * math.comb(-shifted + order - 1, order))
                    assert (second - first) % modulus == 0
                    count += 1
    assert count == 8100
    return {'identities': count, 'degrees': '0 <= d <= 8',
            'moduli': '1 <= M <= 20', 'arguments': '-4 <= m <= 4',
            'basis_orders': '0 <= v <= d', 'step': 'M*lcm(1,...,max(1,d))'}


def rna_determinant_identities(examples):
    values = [1]
    for degree in range(1, 100):
        previous = values[degree - 1]
        preceding = values[degree - 2] if degree >= 2 else 0
        convolution = sum(values[index] * values[degree - 2 - index]
                          for index in range(degree - 1))
        values.append(previous - preceding + convolution)
    count = 0
    for row in examples:
        for argument, expected in enumerate(row['values_m_0_through_5']):
            size = row['residue'] + 8 * argument
            assert hankel(values, 3, size) == expected
            count += 1
    assert count == 48
    return {'identities': count, 'coefficient_series': 'x^2*A^2-(1-x+x^2)*A+1=0; A(0)=1',
            'coverage': 'third coefficient shift; all residues modulo 8; 0 <= m <= 5',
            'identity': 'Delta_j^(4)(Phi_1) = det(a_(3+i+k))_(i,k<j)'}


def check():
    counts = {'cycle_parameters': 0, 'normal_matrices': 0, 'orthogonality_equations': 0,
              'boundary_determinants': 0, 'endpoint_polynomial_sequences': 0}
    for parameter in range(1, 41):
        data = fraction_data(parameter)
        assert len(data) == (3 if parameter == 1 else 6 * parameter - 4)
        assert sum(order + 1 for order, _, _ in data) == 2 * parameter * (parameter + 1)
        matrix = [[1, 0], [0, 1]]
        running_sign, normal_sign = -1, 1
        for order, sign, denominator in data:
            leading = coefficient(denominator, order + 1)
            matrix = [[leading * matrix[0][column] - sign * matrix[1][column] for column in range(2)], matrix[0]]
            running_sign *= sign
            normal_sign *= (-1) ** (order * (order + 1) // 2) * running_sign ** (order + 1)
        assert matrix == [[1, 0], [0, 1]]
        assert running_sign == -1 and normal_sign == (-1) ** parameter
        counts['cycle_parameters'] += 1
    for parameter in range(1, 5):
        for block in [1, 2]:
            normal_size = block * period(parameter)
            first, second, _, _ = continuants(parameter, block)
            values = moments(parameter, parameter + 1 + 2 * normal_size)
            assert hankel(values, parameter + 1, normal_size) == 1
            counts['normal_matrices'] += 1
            for exponent in range(normal_size):
                assert sum(value * values[parameter + 1 + index + exponent] for index, value in enumerate(first)) == 0
                assert sum(value * values[parameter + 1 + index + exponent] for index, value in enumerate(second)) == (-1 if exponent == normal_size - 1 else 0)
                counts['orthogonality_equations'] += 2
        for extra_shift in [2, 3, 4]:
            for size in range(21):
                block = (size + extra_shift + period(parameter) - 1) // period(parameter)
                first, second, _, _ = continuants(parameter, block)
                dimension = block * period(parameter) - size
                values = moments(parameter, parameter + 2 + extra_shift + 2 * size)
                actual = hankel(values, parameter + 1 + extra_shift, size)
                assert actual == boundary(first, second, extra_shift, size, dimension)
                counts['boundary_determinants'] += 1
    for parameter in range(1, 7):
        endpoints = [continuants(parameter, block) for block in range(1, 8)]
        for depth in range(5):
            for endpoint in range(4):
                sequence = [coefficient(values[endpoint], depth) for values in endpoints]
                for _ in range(depth + 1):
                    sequence = differences(sequence)
                assert all(value == 0 for value in sequence)
                counts['endpoint_polynomial_sequences'] += 1
    examples = []
    parameter, shift = 1, 4
    values = moments(parameter, 2 * (6 * period(parameter)) + shift)
    for residue in range(period(parameter)):
        sequence = [hankel(values, shift, residue + block * period(parameter)) for block in range(6)]
        newton = []
        remainder = sequence
        while remainder:
            newton.append(remainder[0])
            remainder = differences(remainder)
        assert all(value == 0 for value in newton[2:])
        examples.append({'residue': residue, 'newton_coefficients': newton[:2], 'values_m_0_through_5': sequence})
    return {'schema_version': 1, 'arithmetic': 'exact integers; fraction-free elimination',
            'coverage': {'cycle_identities': '1 <= n <= 40',
                         'normal_and_orthogonality': '1 <= n <= 4; K = 1, 2',
                         'boundary_identity': '1 <= n <= 4; a = 2, 3, 4; 0 <= j <= 20',
                         'endpoint_laws': '1 <= n <= 6; 0 <= h <= 4; 1 <= K <= 7',
                         'example': 'n = 1, ell = 4; all eight residues; 0 <= m <= 5'},
            'counts': counts, 'golden_shift_4': examples,
            'golden_polynomial_identities': golden_polynomial_identities(),
            'modular_newton_identities': modular_newton_identities(),
            'rna_determinant_identities': rna_determinant_identities(examples),
            'interpretation': 'Finite exact identities in the stated ranges. The 220 golden identities and the written degree bound prove the eight affine laws for all m >= 0. The other all-parameter arguments are written in note.tex.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path)
    arguments = parser.parse_args()
    result = check()
    result['checker_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    text = json.dumps(result, indent=2) + '\n'
    if arguments.output:
        arguments.output.write_text(text)
    else:
        print(text, end='')


if __name__ == '__main__':
    main()
