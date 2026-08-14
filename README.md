## Article:

Using combinatorial chemistry and high-throughput screening, we developed a mass tag, called [PSMtag](https://www.biorxiv.org/content/10.1101/2025.05.22.655509v1), with the potential to allow simultaneously increasing proteome coverage and sample throughput for DIA workflows without compromising quantitative accuracy. We open-sourced the synthesis of reagents and [the software](https://github.com/ParallelSquared/jmod/) to allow rapid strides towards 100-1000x increases in throughput to serve the biological community. 

<h2 style="letter-spacing: 2px; font-size: 26px;" id="data">

Data:

</h2>

All raw and processed data from the [article](https://www.biorxiv.org/content/10.1101/2025.05.22.655509v1) are organized in this data repositories: [MassIVE MSV000097968](http://massive.ucsd.edu/ProteoSAFe/status.jsp?task=4906b9cefe1c4aa29bfca03a27fb6a65) and [ProteomeXchange PXD064191](https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD064191). 

<h2 style="letter-spacing: 2px; font-size: 26px;" id="code">
Code:
</h2>

This repository contains the [Python](https://github.com/ParallelSquared/tag/tree/main/code/Python) and [R](https://github.com/ParallelSquared/tag/tree/main/code/R) scripts used to process data and generate figures for the article.

<details>
<summary>Python instructions</summary>
  <h3>Installing dependencies</h3>
  Begin by installing the required packages:
  
  ```
  pip install -r requirements.txt
  ```

  <h3>Producing DDA comparison plots (sage_dda.py)</h3>
  This script reads Sage search results and renders plots comparing tagged vs. untagged peptide scoring and fragmentation. It takes arguments in the following format:
  
  ```
  python3 sage_dda.py [LF_SAGE_RESULTS_PATH] [TAGGED_SAGE_RESULTS_PATH]
  ```
  where LF_SAGE_RESULTS_PATH and TAGGED_SAGE_RESULTS_PATH are paths to Sage output directories, which must contain:
  - results.sage.tsv
  - lfq.tsv
  - matched_fragments.tsv

To replicate publication figures, use the following Sage output directories provided in the MassIVE FTP repository MSV000097968:
- search/search/LF_28 <i>(Tryptic label-free results, NCE28, Astral)</i>
- search/search/LF_24 <i>(Tryptic label-free results, NCE24, Astral)</i>
- search/search/T6_24 <i>(Tryptic PSMtag results, NCE24, Astral)</i>
   
</details>

<h2 style="letter-spacing: 2px; font-size: 26px;" id="media">

Media:

</h2>

Miscellaneous information, including publicly-available presentations about the [article](https://www.biorxiv.org/content/10.1101/2025.05.22.655509v1) are available through: [parallelsq.org/PSMtags](https://www.parallelsq.org/psmtags).
