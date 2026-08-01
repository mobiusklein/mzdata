import gzip
import io
import itertools
import json
import re
import sys
from enum import IntFlag

import fastobo
from fastobo.doc import OboDoc
from fastobo.id import PrefixedIdent
from fastobo.term import (
    DefClause,
    IsAClause,
    NameClause,
    TermFrame,
)

ACQUISITION_SW = PrefixedIdent("MS", "1001455")
ANALYSIS_SW = PrefixedIdent("MS", "1001456")
DP_SW = PrefixedIdent("MS", "1001457")

segment_pattern = re.compile(r"(_[a-zA-Z])")

class SoftwareType(IntFlag):
    NoType = 0
    Analysis = 0b00000001
    DataProcessing = 0b00000010
    Acquisition = 0b00000100


def collect_software_types(cv: OboDoc) -> tuple[set[PrefixedIdent], dict[PrefixedIdent, TermFrame]]:
    term: TermFrame
    id_to_clause = {}
    software_ids = {
        PrefixedIdent("MS", "1000531")
    }
    for term in itertools.chain(cv, cv):
        id_to_clause[term.id] = term
        for clause in term:
            if isinstance(clause, IsAClause) and clause.term in software_ids:
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
    if "π" in vname:
        vname = vname.replace("π", "pi")

    vname: str = segment_pattern.sub(format_name, vname.replace(" ", "_").replace("software", "Software"))
    vname: str = vname[0].upper() + vname[1:]

    if vname[0].isdigit():
        vname = "_" + vname

    return f"""
    #[term(cv=MS, accession={term.id.local}, name="{name}", flags={{{int(flags)}}}, parents={{{json.dumps(parents)}}})]
    #[doc="{name} - {descr}"]
    {vname},"""


def generate_term_enum(terms: list[TermFrame]):
    buffer = io.StringIO()
    buffer.write("pub enum SoftwareTerm {")
    for term in terms:
        buffer.write(make_entry_for(term))
    buffer.write("\n}")
    return buffer.getvalue()


def main():
    with gzip.open("./cv/psi-ms.obo.gz") as handle:
        cv: OboDoc = fastobo.load(handle)
        software_ids, id_to_clause = collect_software_types(cv)
        sw_terms = list(map(id_to_clause.get, sorted(software_ids)))
        text = generate_term_enum(sw_terms).encode('utf8')
        sys.stdout.buffer.write(text)


if __name__ == "__main__":
    main()