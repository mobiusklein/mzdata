import argparse
import gzip
import json
import io
import itertools
import re

from enum import IntFlag
from typing import Tuple, Dict, Set, List

import fastobo
from fastobo.term import TermFrame, IsAClause, NameClause, RelationshipClause, DefClause

from fastobo.doc import OboDoc

from fastobo.id import PrefixedIdent

segment_pattern = re.compile(r"(_[a-zA-Z])")

class ValueType(IntFlag):
    NoType = 0
    String = 0b00000001
    Integer = 0b00000010
    Float = 0b00000100
    Double = 0b00001000
    NonNegativeInteger = 0b00010000
    PositiveInteger = 0b00100000
    DateTime = 0b01000000
    Boolean = 0b10000000

    ListOf = 0b1000000000000000


xsd_to_type = {
    "xsd:int": ValueType.Integer,
    "xsd:integer": ValueType.Integer,
    "xsd:string": ValueType.String,
    "xsd:float": ValueType.Float,
    "xsd:double": ValueType.Double,
    "xsd:nonNegativeInteger": ValueType.NonNegativeInteger,
    "xsd:positiveInteger": ValueType.PositiveInteger,
    "xsd:dateTime": ValueType.DateTime,
    "xsd:boolean": ValueType.Boolean,
}


COMPONENT_TO_ENUM = {
    "mass-analyzer": 'MassAnalyzer',
    "ionization-type": 'IonizationType',
    "inlet-type": "InletType",
    "detector-type": "DetectorType",
    "collision-energy": "CollisionEnergy"
}


COMPONENT_TO_TERM = {
    "mass-analyzer": PrefixedIdent("MS", "1000443"),
    "ionization-type": PrefixedIdent("MS", "1000008"),
    "inlet-type": PrefixedIdent("MS", "1000007"),
    "detector-type": PrefixedIdent("MS", "1000026"),
    "collision-energy": PrefixedIdent("MS", "1000045"),
}


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "component",
        choices=[
            "mass-analyzer",
            "ionization-type",
            "inlet-type",
            "detector-type",
            "collision-energy",
            "-",
        ],
    )
    parser.add_argument("-c", "--curie")
    parser.add_argument("-t", "--type-name")
    return parser


def collect_components(
    cv: OboDoc,
    base_term: PrefixedIdent
) -> Tuple[Set[PrefixedIdent], Dict[PrefixedIdent, TermFrame]]:
    term: TermFrame
    id_to_clause = {}
    component_ids = {base_term}
    # Make multiple passes
    for term in itertools.chain(cv, cv):
        id_to_clause[term.id] = term
        for clause in term:
            if isinstance(clause, IsAClause):
                if clause.term in component_ids:
                    component_ids.add(term.id)
    return component_ids, id_to_clause


def format_name(match: re.Match) -> str:
    return match.group(1)[-1].upper()


def find_name(term: TermFrame):
    for clause in term:
        if isinstance(clause, NameClause):
            name = str(clause.name)
            return name
    else:
        raise LookupError(f"Term name not found for {term.id!s}")


def make_entry_for(term: TermFrame):
    name = None
    flags = ValueType.NoType
    descr = ""
    parents = []
    for clause in term:
        if isinstance(clause, NameClause):
            name = str(clause.name)
        if isinstance(clause, IsAClause):
            parents.append(str(clause.term))
        if isinstance(clause, RelationshipClause):
            if str(clause.typedef) == 'has_value_type':
                flags |= xsd_to_type[str(clause.term)]
        if isinstance(clause, DefClause):
            descr = re.sub(
                r"(\[|\])",
                lambda m: "\\\\" + m.group(1),
                str(clause.definition).replace('"', "'"),
            )

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

    if vname[0].isdigit():
        vname = "_" + vname

    return f"""
    #[term(cv=MS, accession={term.id.local}, name="{name}", flags={{{int(flags)}}}, parents={{{json.dumps(parents)}}})]
    #[doc="{name} - {descr}"]
    {vname},"""


def generate_term_enum(terms: List[TermFrame], type_name: str):
    buffer = io.StringIO()
    buffer.write("pub enum $Term {".replace("$", type_name))
    for term in terms:
        buffer.write(make_entry_for(term))
    buffer.write("\n}")
    return buffer.getvalue()


def main():
    parser = make_parser()
    args = parser.parse_args()
    component = args.component

    if component == '-':
        term = PrefixedIdent(*args.curie.split(":"))
        type_name = args.type_name
    else:
        term = COMPONENT_TO_TERM[component]
        type_name = COMPONENT_TO_ENUM[component]

    cv: OboDoc = fastobo.load(gzip.open("./cv/psi-ms.obo.gz"))
    component_ids, id_to_clause = collect_components(cv, term)
    if type_name is None:
        t = find_name(id_to_clause[term])
        type_name = t.title().replace(" ", "")

    term_specs = list(map(id_to_clause.get, sorted(component_ids)))
    text = generate_term_enum(term_specs, type_name)
    print(text)


if __name__ == "__main__":
    main()