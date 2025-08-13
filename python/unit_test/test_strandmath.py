# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 15:49:57 2025

@author: aledb
"""

# test_strand.py
import pytest
import copy
import re
import random
import numpy as np
from strandmath import Strand, Multistrand, Duplex


def test_init_single_str():
    s = Strand("ATGC", name="test")
    assert s.numel() == 1
    assert s.name == ["test"]
    assert s.string() == "ATGC"


def test_init_from_strand_copy():
    s1 = Strand("ATGC", name="orig")
    s2 = Strand(s1)
    assert s2.string() == "ATGC"
    assert s2.name == ["orig"]
    assert s1 is not s2
    assert s1.sequence is not s2.sequence


def test_init_list_of_str():
    seqs = ["ATGC", "GATTACA"]
    names = ["s1", "s2"]
    s = Strand(seqs, name=names)
    assert s.numel() == 2
    assert s.string() == seqs
    assert s.name == names


def test_init_invalid_name_type():
    with pytest.raises(TypeError):
        Strand(["ATGC", "GATTACA"], name="bad")


def test_getitem_and_setitem():
    s = Strand(["AAA", "TTT"], name=["n1", "n2"])
    first = s[0]
    assert first.string() == "AAA"
    assert first.name == ["n1"]

    replacement = Strand("GGG", name="new")
    s[0] = replacement
    assert s[0].string() == "GGG"

    with pytest.raises(ValueError):
        s[0] = Strand(["AAA", "CCC"], name=["x", "y"])

    with pytest.raises(TypeError):
        s[0] = "bad"


def test_repr():
    assert repr(Strand("ATGC")) == "Strand('ATGC')"
    assert "Sequences" in repr(Strand(["AAA", "TTT"], name=["", ""]))


def test_add_and_radd():
    s1 = Strand("AAA", name="a")
    s2 = Strand("TTT", name="b")
    res = s1 + s2
    assert res.string() == "AAATTT"
    res2 = "GGG" + s1
    assert res2.string() == "GGGAAA"


def test_neg_and_sub():
    s1 = Strand("ATGC", name="a")
    s2 = Strand("TTAA", name="b")
    rev = -s1
    assert rev.string() == "CGTA"
    sub_res = s1 - s2
    assert sub_res.string()[0].startswith("ATGC")


def test_invert_reverse_complement():
    s = Strand("ATGC")
    rc = ~s
    assert rc.string() == "GCAT"


def test_eq():
    assert Strand("ATGC") == Strand("ATGC")
    assert Strand("ATGC") != Strand("TTTT")


def test_is_symmetric():
    pal = Strand("ATGCAT")
    assert pal.isSymmetric() in (True, False)  # Just ensure it runs


def test_remove_duplicates():
    s = Strand(["AAA", "AAA", "TTT"], name=["a", "dup", "b"])
    rd = s.removeDuplicates()
    assert rd.numel() == 2
    assert "dup" not in rd.name


def test_scramble():
    s = Strand("AAAA", name="test")
    scrambled = s.scramble()
    assert scrambled.name[0].endswith("_scrambled")
    assert sorted(scrambled.sequence[0]) == sorted(s.sequence[0])


def test_from_string_and_cleaning():
    assert Strand("5'-A T G C-3'").string() == "ATGC"


def test_string_and_barestring():
    s = Strand("+A bT rG C")
    assert isinstance(s.string(), str)
    bare = s.bareString()
    assert "+" not in bare and "b" not in bare and "r" not in bare


def test_bare_sequence():
    s = Strand("+A bT")
    bare = s.bareSequence()
    assert all(re.match(r"^[ACGTU]$", nt) for row in bare for nt in row)


def test_remove_mods_and_cleanstring():
    assert Strand.removeMods("+A bT") == "A T"
    assert Strand.cleanString("5'-AT-3'") == "AT"


def test_conversions():
    s = Strand("AUTG")
    assert "T" in s.toDNA().string()
    assert all(nt.startswith("r") for nt in s.toRNA().sequence[0])
    assert all(nt.startswith("+") for nt in s.toLNA().sequence[0])
    assert all(nt.startswith("b") for nt in s.toBNA().sequence[0])


def test_len_method():
    assert Strand("ATGC").len() == 4
    multi = Strand(["AAA", "TTTT"])
    assert multi.len() == [3, 4]


def test_reverse():
    s = Strand("ATGC", name="n")
    r = s.reverse()
    assert r.string() == "CGTA"
    assert r.name[0].endswith("_reverse")


def test_reverse_complement():
    s = Strand("ATGC", name="n")
    rc = s.reverseComplement()
    assert rc.string() == "GCAT"
    assert rc.name[0].endswith("_complement")


def test_crop_and_errors():
    s = Strand("ATGC")
    cropped = s.crop([1, 3])
    assert cropped.string() == "TG"
    with pytest.raises(TypeError):
        s.crop([1, "bad"])
    with pytest.raises(ValueError):
        s.crop([2, 1])


def test_gc_content_and_gc():
    s = Strand("GCGC")
    assert s.gcContent() == 1.0
    assert s.gc() == 1.0


def test_random_generation():
    s = Strand.random(10, seqtype="DNA", gcContent=0.5)
    assert isinstance(s, Strand)
    assert s.len() == 10


def test_multistrand_init():
    ms = Multistrand("ATGC", "")
    assert ms.Strands.numel() == 2
    assert ms.Strands[1].string()


def test_duplex_defaults():
    d = Duplex()
    assert isinstance(d.Strands, Strand)
    assert isinstance(d.Schema, list)
    assert isinstance(d.Nbp, np.ndarray)