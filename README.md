# MTase-classification
 Code for detection MTase class and regions of subdomains and motifs

For running scripts look at `Classification-pipeline.ipynb`.

## Classification algorithm

Classification of MTases performed by scripts `classification_SUPFAM.py` and `Classification_Pfam.py`. 

For cat-profile from SUPFAM the algorithm assigns a class according to the following principle:


<table>
  <tr>
   <td><strong>Class</strong>
   </td>
   <td><strong>Criteria for assigning a catalytic domain to a class</strong>
   </td>
  </tr>
  <tr>
   <td>A
   </td>
   <td>1. Profiles from group A find cat-motif (4 letters found) and sam-motif in sequence (8 letters found);
<p>
2. The cat-motif is characteristic of Class A,  or has not been encountered before;
<p>
3. In the МТase sequence, first comes the sam- motif, then cat-motif;
<p>
4. The profile finding includes at least 2 profile positions corresponding to the beginning of the SAM-binding subdomain (i.e. the Hd1 helix).
   </td>
  </tr>
  <tr>
   <td>B
   </td>
   <td>1. Profiles from group B find cat-motif (4 letters found) and sam-motif in sequence (8 letters found);
<p>
2. The cat-motif is characteristic of Class B,  or has not been encountered before;
<p>
3. In the МТase sequence, first comes the cat- motif, then sam-motif;
<p>
4. The region of the profile S7-Hu3 is aligned at least 10 positions. The region of the profile S4-Hu2-S3-Hu1-S2-S1 is aligned at least 20 positions.  
<p>
5. The sequence of profile hits corresponds to the sequence S7-Hu3-S4-Hu2-S3-Hu1-S2-S1-TRD-Hd1-S5-Hd2-S6-Hd3.
   </td>
  </tr>
  <tr>
   <td>D
   </td>
   <td>1. Profiles from group D find cat-motif (4 letters found) and sam-motif in sequence (8 letters found);
<p>
2. The cat-motif is characteristic of Class D,  or has not been encountered before;
<p>
3.  In the МТase sequence, first comes the sam- motif, then cat-motif;
<p>
4. The region of the profile S7-Hu3 is aligned at least 10 positions. The region of the profile Hd1-S5-Hd2-S6-Hd3  profile is aligned at least 20 positions;
<p>
5. The sequence of profile hits corresponds to the sequence Hd1-S5-Hd2-S6-Hd3-TRD-S7-Hu3-S4-Hu2-S3-Hu1-S2-S1.
   </td>
  </tr>
  <tr>
   <td>C
   </td>
   <td>1. Profiles from group C find cat-motif (4 letters found) and sam-motif in sequence (8 letters found);
<p>
2. The cat-motif is characteristic of Class C,  or has not been encountered before;
<p>
3. In the МТase sequence, first comes the sam- motif, then cat-motif;
<p>
4. The region of the profile Hd1 is aligned at least 10 positions. The region of the profile S5-Hd2-S6-Hd3-S7-Hu3  profile is aligned at least 20 positions;
<p>
5. The sequence of profile hits corresponds to the sequence S5-Hd2-S6-Hd3-S7-Hu3-S4-Hu2-S3-Hu1-S2-S1-TRD-Hd1.
   </td>
  </tr>
  <tr>
   <td>K
   </td>
   <td>1. Profiles from group C find cat-motif (4 letters found) and sam-motif in sequence (8 letters found);
<p>
2. The cat-motif is characteristic of Class C,  or has not been encountered before;
<p>
3. In the МТase sequence, first comes the sam- motif, then cat-motif;
<p>
4. The sequence of profile hits corresponds to the sequence Hd1-S5-Hd2-S6-Hd3-S7-Hu3-S4-Hu2-S3-Hu1-S2-S1;
   </td>
  </tr>
  <tr>
   <td>L
   </td>
   <td>1 Profiles from group B find cat-motif (4 letters found) and sam-motif in sequence (8 letters found);
<p>
2 The cat-motif is characteristic of class B,  or has not been encountered before;
<p>
3  In the МТase sequence, first comes the sam- motif, then cat-motif;
<p>
4. The sequence of profile hits corresponds to the sequence Hd1-S5-Hd2-S6-Hd3-S7-Hu3-S4-Hu2-S3-Hu1-S2-S1;
<p>
5 The distance between regions Hd1-S5-Hd2-S6-Hd3  and S7-Hu3  does not exceed 20 a.c.
   </td>
  </tr>
  <tr>
   <td>M
   </td>
   <td>1 Profiles from group B find cat-motif (4 letters found) and sam-motif in sequence (8 letters found);
<p>
2 The cat-motif is characteristic of class B,  or has not been encountered before;
<p>
3  In the МТase sequence, first comes the sam- motif, then cat-motif;
<p>
4. The sequence of profile hits corresponds to the sequence Hd1-S5-Hd2-S6-Hd3-S7-Hu3-S4-Hu2-S3-Hu1-S2-S1;
<p>
5 The distance between regions Hd1-S5-Hd2-S6-Hd3  and S7-Hu3  does not exceed 20 a.c.../
   </td>
  </tr>
  <tr>
   <td>I
   </td>
   <td>1. Profiles from group B find cat-motif (4 letters found);
<p>
2. Expert assignment.
   </td>
  </tr>
</table>


For cat-profile from PFAM the algorithm assigns a class according to the following principle:


<table>
  <tr>
   <td><strong>Class</strong>
   </td>
   <td><strong>Criteria for assigning a catalytic domain to a class</strong>
   </td>
  </tr>
  <tr>
   <td>E
   </td>
   <td>1. Profiles from group E find cat-motif (4 letters found) and sam-motif in sequence (more than 5 letters found);
<p>
2. The cat-motif is NPPY,  or has not been encountered before;
<p>
3. In the МТase sequence, first comes the sam- motif, then cat-motif.
   </td>
  </tr>
  <tr>
   <td>F
   </td>
   <td>1. Profiles from group F find cat-motif (4 letters found) and sam-motif in sequence (2 letters found);
<p>
2. The cat-motif is NPPF,  or has not been encountered before;
<p>
3. In the МТase sequence, first comes the sam- motif, then cat-motif.
   </td>
  </tr>
  <tr>
   <td>G
   </td>
   <td>1. Profiles from group F find cat-motif (4 letters found) and sam-motif in sequence (8 letters found);
<p>
2. The cat-motif is DPPW,  or EPPW or  has not been encountered before;
<p>
3. In the МТase sequence, first comes the cat- motif, then sam-motif. 
   </td>
  </tr>
</table>

