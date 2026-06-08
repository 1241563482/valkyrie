import numpy as np


qpoint = 16
qband = 3
nsig = 30


def _extract_columns():
    dos_rows = []
    lambda_values = []
    with open(f"elph.inp_lambda.{qpoint}") as file:
        for lineno, line in enumerate(file, start=1):
            parts = line.split()
            if "Ef" in line:
                dos_rows.append([lineno, parts[2], parts[7]])
            if "lambda" in line:
                lambda_values.append(parts[2])
    np.savetxt("DOS_Ef.txt", np.array(dos_rows, dtype=object), fmt="%s")
    np.savetxt("tb.txt", np.array(lambda_values, dtype=object), fmt="%s")


def main():
    _extract_columns()
    qk_lambda = np.loadtxt('tb.txt')

    sumlambda = []
    for ns in range(nsig):
        s = 0
        for qb in range(qband):
            s += qk_lambda[ns * qband + qb]
        sumlambda.append(s)
    np.savetxt('lambda_vs_nsig' + str(qpoint) + '.dat', sumlambda)
    print(sumlambda)


if __name__ == "__main__":
    main()
