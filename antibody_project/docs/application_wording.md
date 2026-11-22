# Wording for the application

## Filled in (pipeline run 2026-11-22) — Servier Symphogen, Predictive Antibody Purification

CV line:

> **Sequence-based purification descriptors for clinical-stage antibodies** (Python, open
> source; `github.com/bozmen17/antibody-purification-descriptors`) — Built a pipeline computing pI, pH-dependent charge,
> charge/hydrophobic-patch proxies and CDR liability motifs for 1,164 Fv arms of 1,069
> therapeutics in Thera-SAbDab. PCA and clustering showed purification-relevant charge
> variation is spread across all formats rather than following them (ARI ≈ 0.02); a
> bispecific module quantifies arm-to-arm charge separation from homodimer and
> light-chain mispairs, flagging ~1 in 4 IgG-like bispecifics as having little
> Fv-level charge handle against their homodimers.

Cover-letter sentence:

> To get closer to your problem before applying, I built an open-source pipeline computing
> purification-relevant sequence descriptors for the 1,069 clinical-stage antibodies in
> Thera-SAbDab, including an arm-to-arm charge calculation that estimates, from sequence
> alone, how much separation a bispecific offers against its homodimer and light-chain
> mispairs — about a quarter of IgG-like bispecifics have less than one elementary charge
> of Fv-level handle — and connecting such predictions to measured chromatographic
> behaviour is exactly the problem this position addresses.

---

Placeholders in `<angle brackets>`. Fill the numbers in only after you have run the
pipeline — a figure you can't defend in an interview is worse than no figure.

## CV line

Under a "Projects" heading, or at the end of the relevant experience block:

> **In silico antibody developability descriptors** — Python, open source.
> Built a sequence-only descriptor pipeline over the <n> clinical-stage antibody
> therapeutics in Thera-SAbDab (isoelectric point and charge curves at process pH,
> charge- and hydrophobic-patch proxies, CDR liability motifs), then used PCA and
> clustering to test how far purification-relevant properties separate antibody
> formats. `github.com/bozmen17/antibody-purification-descriptors`

Shorter, if space is tight:

> **In silico antibody developability descriptors** (Python; `github.com/bozmen17/antibody-purification-descriptors`) —
> charge, patch and liability descriptors for <n> clinical-stage antibodies from
> Thera-SAbDab, with PCA and clustering by molecular format.

## Cover-letter sentence

Pick one; both do the same job of moving you from "I studied this" to "I have started
on your problem".

> To get closer to the problem before applying, I built a small open-source pipeline
> that computes purification-relevant descriptors — isoelectric point, charge at
> process pH, and charge- and hydrophobic-patch proxies — for every clinical-stage
> antibody in Thera-SAbDab, and used PCA and clustering to test how far these
> properties track molecular format; the natural next step, and what attracts me to
> this role, is to connect such descriptors to measured chromatographic behaviour.

> I have already started on this class of problem: a small open-source project of mine
> computes charge and hydrophobicity descriptors for the clinical-stage antibody
> landscape from public sequence data and examines their structure by PCA and
> clustering, which left me wanting exactly what this position offers — in-house
> process data to turn descriptors into predictions.

## If you get a technical interview, be ready for

- **"Why those descriptors?"** Because the separations in a platform process resolve on
  charge and on exposed apolar surface, so pI, the charge curve across the pH range the
  process uses, and patch proxies are the sequence-level quantities closest to what the
  column sees.
- **"Why is the agreement with format weak?"** Thera-SAbDab gives variable domains;
  format is mostly a constant-region, glycan and linker decision. Know this answer
  before the question comes — it is the main limitation of the study and you want to be
  the one who raises it.
- **"What would you do with real data?"** Supervised model on a measured endpoint —
  HIC retention time, Protein A elution pH, or yield at a defined step — with the
  descriptors as features; and structure-based patches on homology models instead of
  windowed sums.
- **"Why windowed sums rather than real patches?"** Honest answer: sequence is all the
  database gives for every entry, and a window is the cheapest stand-in. Say what the
  correct version is (structure- or model-based patch calculation).

## Remaining time, if you want to go further

With two months, in order of return on effort:

1. Run it, read the output, and write the Results section of the README. Non-negotiable.
2. Swap the heuristic loop boundaries for ANARCI numbering. Removes the weakest part
   of the method and is a day's work.
3. One supervised result on a published dataset with a measured endpoint. This is what
   turns the project from descriptive into predictive, and it is the single thing most
   likely to be discussed in the interview.
4. Only then: polish, a notebook walkthrough, or a short write-up on your website.
