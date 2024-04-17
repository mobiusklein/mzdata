import argparse
import gzip
import json
import io
import re

from enum import IntFlag
from typing import Tuple, Dict, Set, List

import fastobo
from fastobo.term import (
    TermFrame,
    IsAClause,
    NameClause,
)

from fastobo.doc import OboDoc

from fastobo.id import PrefixedIdent

segment_pattern = re.compile(r"(_[a-zA-Z])")

MASS_ANALYZER_BASE_TERM = PrefixedIdent("MS", "1000443")
IONIZATION_BASE_TERM = PrefixedIdent("MS", "1000008")


COMPONENT_TO_ENUM = {
    "mass-analyzer": 'MassAnalyzer',
    "ionization-type": 'IonizationType',
    "inlet-type": "InletType",
    "detector-type": "DetectorType",
}


COMPONENT_TO_TERM = {
    "mass-analyzer": MASS_ANALYZER_BASE_TERM,
    "ionization-type": IONIZATION_BASE_TERM,
    "inlet-type": PrefixedIdent("MS", "1000007"),
    "detector-type": PrefixedIdent("MS", "1000026"),
}


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("component", choices=["mass-analyzer", "ionization-type", "inlet-type", "detector-type"])
    return parser


def collect_components(
    cv: OboDoc,
    component: str,
) -> Tuple[Set[PrefixedIdent], Dict[PrefixedIdent, TermFrame]]:
    base_term = COMPONENT_TO_TERM[component]
    term: TermFrame
    id_to_clause = {}
    component_ids = {base_term}
    for term in cv:
        id_to_clause[term.id] = term
        for clause in term:
            if isinstance(clause, IsAClause):
                if clause.term in component_ids:
                    component_ids.add(term.id)
    return component_ids, id_to_clause


def format_name(match: re.Match) -> str:
    return match.group(1)[-1].upper()


def make_entry_for(term: TermFrame):
    name = None
    flags = 0
    parents = []
    for clause in term:
        if isinstance(clause, NameClause):
            name = str(clause.name)
        if isinstance(clause, IsAClause):
            parents.append(str(clause.term))

    vname = name
    if "-" in vname:
        vname = vname.replace("-", "_")
    if ":" in vname:
        vname = vname.replace(":", "_")
    if "/" in vname:
        vname = vname.replace("/", "_")
    if "+" in vname:
        vname = vname.replace("+", "plus")
    if "!" in vname:
        vname = vname.replace("!", "_")

    vname = segment_pattern.sub(
        format_name, vname.replace(" ", "_")
    )
    vname = vname[0].upper() + vname[1:]

    return f"""
    #[term(cv=MS, accession={term.id.local}, name="{name}", flags={{{int(flags)}}}, parents={{{json.dumps(parents)}}})]
    {vname},"""


def generate_term_enum(terms: List[TermFrame], component: str):
    buffer = io.StringIO()
    buffer.write("pub enum $Term {".replace("$", COMPONENT_TO_ENUM[component]))
    for term in terms:
        buffer.write(make_entry_for(term))
    buffer.write("\n}")
    return buffer.getvalue()


def main():
    parser = make_parser()
    args = parser.parse_args()
    component = args.component
    cv: OboDoc = fastobo.load(gzip.open("./cv/psi-ms.obo.gz"))
    component_ids, id_to_clause = collect_components(cv, component)
    term_specs = list(map(id_to_clause.get, sorted(component_ids)))
    text = generate_term_enum(term_specs, component)
    print(text)


if __name__ == "__main__":
    main()