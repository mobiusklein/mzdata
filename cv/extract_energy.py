import gzip
import json
import io
import itertools
import re

import sys
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

ROOT_TERM = PrefixedIdent("MS", "1000045")
EXTRA_ROOTS = [
    PrefixedIdent("MS", "1000138"),
    PrefixedIdent("MS", "1002680"),
    PrefixedIdent("MS", "1003410")
]

segment_pattern = re.compile(r"(_[a-zA-Z])")


def collect_components(
    cv: OboDoc, base_term: PrefixedIdent
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
    parents = []
    descr = ""
    for clause in term:
        if isinstance(clause, NameClause):
            name = str(clause.name)
        if isinstance(clause, IsAClause):
            parents.append(str(clause.term))
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

    vname = segment_pattern.sub(format_name, vname.replace(" ", "_"))
    vname = vname[0].upper() + vname[1:]

    if vname[0].isdigit():
        vname = "_" + vname

    return f"""
    #[term(cv=MS, accession={term.id.local}, name="{name}", flags={{{0}}}, parents={{{json.dumps(parents)}}})]
    #[doc = "{name} - {descr}"]
    {vname}(f32),"""


def generate_term_enum(terms: List[TermFrame], type_name: str):
    buffer = io.StringIO()
    buffer.write("pub enum $Term {".replace("$", type_name))
    for term in terms:
        buffer.write(make_entry_for(term))
    buffer.write("\n}")
    return buffer.getvalue()


def merge_term_sets(term_sets: List[Tuple[Set, Dict]]) -> Tuple[Set, Dict]:
    base_term_ids, base_id_to_clause = map(lambda x: x.copy(), term_sets[0])
    for (term_ids, id_to_clause) in term_sets[1:]:
        base_term_ids.update(term_ids)
        base_id_to_clause.update(id_to_clause)
    return (base_term_ids, base_id_to_clause)


def main():
    cv: OboDoc = fastobo.load(gzip.open("./cv/psi-ms.obo.gz"))
    term_ids, id_to_clause = merge_term_sets([collect_components(cv, root) for root in [ROOT_TERM] + EXTRA_ROOTS])
    # t = find_name(id_to_clause[ROOT_TERM])
    # type_name = t.title().replace(" ", "")
    type_name = "DissociationEnergy"

    term_specs = list(map(id_to_clause.get, sorted(term_ids)))
    text = generate_term_enum(term_specs, type_name).encode('utf8')
    sys.stdout.buffer.write(text)


if __name__ == "__main__":
    main()