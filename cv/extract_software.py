import gzip
import json
import io
import itertools
import re

from enum import IntFlag
from typing import Tuple, Dict, Set, List

import fastobo
from fastobo.term import (
    TermFrame,
    IsAClause,
    NameClause,
    DefClause,
)

from fastobo.doc import OboDoc

from fastobo.id import PrefixedIdent

ACQUISITION_SW = PrefixedIdent("MS", "1001455")
ANALYSIS_SW = PrefixedIdent("MS", "1001456")
DP_SW = PrefixedIdent("MS", "1001457")

segment_pattern = re.compile(r"(_[a-zA-Z])")

class SoftwareType(IntFlag):
    NoType = 0
    Analysis = 0b00000001
    DataProcessing = 0b00000010
    Acquisition = 0b00000100


def collect_software_types(cv: OboDoc) -> Tuple[Set[PrefixedIdent], Dict[PrefixedIdent, TermFrame]]:
    term: TermFrame
    id_to_clause = {}
    software_ids = {
        PrefixedIdent("MS", "1000531")
    }
    for term in itertools.chain(cv, cv):
        id_to_clause[term.id] = term
        for clause in term:
            if isinstance(clause, IsAClause):
                if clause.term in software_ids:
                    software_ids.add(term.id)
    return software_ids, id_to_clause

def format_name(match: re.Match) -> str:
    return match.group(1)[-1].upper()

def make_entry_for(term: TermFrame):
    name = None
    flags = SoftwareType.NoType
    parents = []
    descr = ''
    for clause in term:
        if isinstance(clause, NameClause):
            name = str(clause.name)
        if isinstance(clause, IsAClause):
            parents.append(str(clause.term))
            if clause.term == DP_SW:
                flags |= SoftwareType.DataProcessing
            elif clause.term == ANALYSIS_SW:
                flags |= SoftwareType.Analysis
            elif clause.term == ACQUISITION_SW:
                flags |= SoftwareType.Acquisition
        if isinstance(clause, DefClause):
            descr = re.sub(
                r"(\[|\])",
                lambda m: "\\\\" + m.group(1),
                str(clause.definition).replace('"', "'"),
            )

    vname: str = name
    if "-" in vname:
        vname = vname.replace("-", "_")
    if ":" in vname:
        vname = vname.replace(":", "_")
    if '/' in vname:
        vname = vname.replace('/', '_')
    if "+" in vname:
        vname = vname.replace("+", "plus")
    if "!" in vname:
        vname = vname.replace("!", "_")

    vname: str = segment_pattern.sub(format_name, vname.replace(" ", "_").replace("software", "Software"))
    vname: str = vname[0].upper() + vname[1:]

    if vname[0].isdigit():
        vname = "_" + vname

    return f"""
    #[term(cv=MS, accession={term.id.local}, name="{name}", flags={{{int(flags)}}}, parents={{{json.dumps(parents)}}})]
    #[doc="{name} - {descr}"]
    {vname},"""


def generate_term_enum(terms: List[TermFrame]):
    buffer = io.StringIO()
    buffer.write("pub enum SoftwareTerm {")
    for term in terms:
        buffer.write(make_entry_for(term))
    buffer.write("\n}")
    return buffer.getvalue()


def main():
    cv: OboDoc = fastobo.load(gzip.open("./cv/psi-ms.obo.gz"))
    software_ids, id_to_clause = collect_software_types(cv)
    sw_terms = list(map(id_to_clause.get, sorted(software_ids)))
    text = generate_term_enum(sw_terms)
    print(text)


if __name__ == "__main__":
    main()