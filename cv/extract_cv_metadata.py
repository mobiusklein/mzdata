import gzip
import argparse

import fastobo

from fastobo.doc import OboDoc
from fastobo.header import DataVersionClause



def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('clause', choices=["data-version"])
    return parser.parse_args()


def main():
    args = make_parser()
    cv: OboDoc = fastobo.load(gzip.open("./cv/psi-ms.obo.gz"))

    clause = args.clause

    target_type = None.__class__
    if clause == 'data-version':
        target_type = DataVersionClause

    for clause in cv.header:
        if isinstance(clause, target_type):
            print(clause.raw_value())


if __name__ == "__main__":
    main()