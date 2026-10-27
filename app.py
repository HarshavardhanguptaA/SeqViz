"""
DNA Sequence Alignment Tool
Educational implementation of LCS, Levenshtein edit distance, and Smith-Waterman.
Run: streamlit run app.py
"""
from __future__ import annotations

import html
import time

import pandas as pd
import streamlit as st

from seqviz.dna_algorithms import (
    edit_distance,
    lcs,
    naive_edit_distance,
    smith_waterman,
    validate_dna,
)


def matrix_frame(matrix: list[list[int]], a: str, b: str) -> pd.DataFrame:
    return pd.DataFrame(matrix, index=["∅"] + list(a), columns=["∅"] + list(b))


def alignment_html(a: str, b: str) -> str:
    top, middle, bottom = [], [], []
    for x, y in zip(a, b):
        if x == y:
            top.append(f'<span style="background:#d1fae5;color:#065f46;padding:2px 3px;border-radius:3px">{html.escape(x)}</span>')
            bottom.append(f'<span style="background:#d1fae5;color:#065f46;padding:2px 3px;border-radius:3px">{html.escape(y)}</span>')
            middle.append("|")
        elif x == "-" or y == "-":
            top.append(f'<span style="background:#ffedd5;color:#9a3412;padding:2px 3px;border-radius:3px">{html.escape(x)}</span>')
            bottom.append(f'<span style="background:#ffedd5;color:#9a3412;padding:2px 3px;border-radius:3px">{html.escape(y)}</span>')
            middle.append(" ")
        else:
            top.append(f'<span style="background:#fee2e2;color:#991b1b;padding:2px 3px;border-radius:3px">{html.escape(x)}</span>')
            bottom.append(f'<span style="background:#fee2e2;color:#991b1b;padding:2px 3px;border-radius:3px">{html.escape(y)}</span>')
            middle.append("×")
    return "<div style='font-family:monospace;line-height:2.2;overflow-wrap:anywhere'>" + "".join(top) + "<br>" + "".join(middle) + "<br>" + "".join(bottom) + "</div>"


st.set_page_config(page_title="DNA Sequence Alignment Tool", page_icon="🧬", layout="wide")
st.title("DNA Sequence Alignment Tool")
st.caption("Dynamic programming for DNA comparison · LCS · Levenshtein edit distance · Smith–Waterman")
with st.sidebar:
    st.header("Sequences")
    seq1_raw = st.text_area("DNA sequence 1", "ACGTACGT", height=100)
    seq2_raw = st.text_area("DNA sequence 2", "ACGTCGT", height=100)
    st.caption("Whitespace is ignored. Accepted bases: A, C, G, T, N.")
    st.divider()
    st.subheader("Smith–Waterman scoring")
    match_score = st.number_input("Match reward", value=2, min_value=1, max_value=10)
    mismatch_score = st.number_input("Mismatch penalty", value=-1, min_value=-10, max_value=0)
    gap_score = st.number_input("Gap penalty", value=-2, min_value=-10, max_value=0)

try:
    seq1 = validate_dna(seq1_raw, "Sequence 1")
    seq2 = validate_dna(seq2_raw, "Sequence 2")
except ValueError as exc:
    st.error(str(exc))
    st.stop()

tab1, tab2, tab3, tab4 = st.tabs(["Overview", "Alignment visualizer", "Mutation simulator", "Performance benchmark"])

with tab1:
    st.subheader("Sequence comparison")
    c1, c2, c3, c4 = st.columns(4)
    distance, edp, operations = edit_distance(seq1, seq2)
    common, common_len, lcsp = lcs(seq1, seq2)
    similarity = max(0.0, 100.0 * (1 - distance / max(len(seq1), len(seq2), 1)))
    c1.metric("Edit distance", distance)
    c2.metric("Similarity", f"{similarity:.2f}%")
    c3.metric("LCS length", common_len)
    c4.metric("Sequence lengths", f"{len(seq1)} / {len(seq2)}")
    st.markdown("**Longest common subsequence (LCS)**")
    st.code(common or "(no common characters)")
    st.markdown("**Edit operations (minimum-cost path)**")
    st.dataframe(pd.DataFrame(operations, columns=["Operation", "Sequence 1", "Sequence 2"]), use_container_width=True, hide_index=True)
    st.info("Similarity is defined here as 100 × (1 − Levenshtein distance / maximum sequence length). This is a normalized edit-distance score, not a clinical or evolutionary identity measure.")
    with st.expander("View dynamic-programming matrices"):
        st.markdown("**Levenshtein edit-distance matrix**")
        st.dataframe(matrix_frame(edp, seq1, seq2), use_container_width=True)
        st.markdown("**LCS matrix**")
        st.dataframe(matrix_frame(lcsp, seq1, seq2), use_container_width=True)

with tab2:
    st.subheader("Smith–Waterman local alignment")
    aligned1, aligned2, sw_matrix, sw_score = smith_waterman(seq1, seq2, int(match_score), int(mismatch_score), int(gap_score))
    st.metric("Best local alignment score", sw_score)
    if aligned1:
        st.markdown(alignment_html(aligned1, aligned2), unsafe_allow_html=True)
        matches = sum(x == y for x, y in zip(aligned1, aligned2))
        st.caption(f"Aligned columns: {len(aligned1)} · Exact matches: {matches} · Match fraction: {100 * matches / max(1, len(aligned1)):.2f}%")
    else:
        st.write("No positive-scoring local alignment was found with these scoring parameters.")
    with st.expander("View Smith–Waterman scoring matrix"):
        st.dataframe(matrix_frame(sw_matrix, seq1, seq2), use_container_width=True)
    st.caption("Green = match · Red = substitution · Orange = gap")

with tab3:
    st.subheader("Explore how mutations affect edit distance")
    mutation_type = st.selectbox("Mutation type", ["Substitution", "Insertion", "Deletion"])
    position = st.number_input("Position (0-based for insertion; 1-based for substitution/deletion)", min_value=0, max_value=max(len(seq1), 1), value=1, step=1)
    base = st.selectbox("Inserted/replacement base", ["A", "C", "G", "T", "N"])
    if st.button("Apply mutation to Sequence 1", type="primary"):
        mutated = seq1
        if mutation_type == "Substitution":
            if not seq1:
                st.error("Sequence 1 is empty.")
            else:
                idx = min(max(int(position) - 1, 0), len(seq1) - 1)
                mutated = seq1[:idx] + base + seq1[idx + 1:]
        elif mutation_type == "Insertion":
            idx = min(int(position), len(seq1))
            mutated = seq1[:idx] + base + seq1[idx:]
        else:
            if not seq1:
                st.error("Sequence 1 is empty.")
            else:
                idx = min(max(int(position) - 1, 0), len(seq1) - 1)
                mutated = seq1[:idx] + seq1[idx + 1:]
        before = edit_distance(seq1, seq2)[0]
        after = edit_distance(mutated, seq2)[0]
        st.write(f"**Original:** `{seq1}`")
        st.write(f"**Mutated:** `{mutated}`")
        st.write(f"**Compared with Sequence 2:** edit distance {before} → {after} (change: {after - before:+d})")
        st.write(f"**New normalized similarity:** {max(0.0, 100 * (1 - after / max(len(mutated), len(seq2), 1))):.2f}%")
        st.caption("This mutation is applied to Sequence 1 for this calculation; it does not overwrite your input.")

with tab4:
    st.subheader("Dynamic programming vs naive recursion")
    st.write("Both methods calculate Levenshtein distance. Dynamic programming reuses subproblem results; naive recursion repeatedly recomputes them.")
    max_len = st.slider("Benchmark sequence length (naive recursion grows exponentially)", min_value=2, max_value=12, value=min(8, max(2, min(len(seq1), len(seq2)))))
    bench_a = (seq1 + "ACGT" * 12)[:max_len]
    bench_b = (seq2 + "TGCA" * 12)[:max_len]
    t0 = time.perf_counter()
    dp_dist = edit_distance(bench_a, bench_b)[0]
    dp_time = time.perf_counter() - t0
    t0 = time.perf_counter()
    recursive_dist = naive_edit_distance(bench_a, bench_b)
    recursive_time = time.perf_counter() - t0
    assert dp_dist == recursive_dist, "The two implementations should return the same distance."
    st.dataframe(pd.DataFrame([
        {"Method": "Dynamic programming", "Edit distance": dp_dist, "Runtime (ms)": dp_time * 1000},
        {"Method": "Naive recursion", "Edit distance": recursive_dist, "Runtime (ms)": recursive_time * 1000},
    ]), use_container_width=True, hide_index=True)
    chart = pd.DataFrame({"Method": ["Dynamic programming", "Naive recursion"], "Runtime (ms)": [dp_time * 1000, recursive_time * 1000]})
    st.bar_chart(chart, x="Method", y="Runtime (ms)")
    st.caption("Runtime varies by machine and Python environment. This small-input benchmark is illustrative, not a rigorous performance study.")
