# Superpixel-based-non-parametric-scene-parsing

Scene parsing refers to the task of labeling every pixel in an image with the class label it belongs to. In this paper, we propose a novel scalable non-parametric scene parsing system based on superpixels correspondence. The non-parametric approach requires almost no training and can scale up to datasets with thousands of labels. This involves retrieving a set of images similar to the query image, followed by superpixel matching of the query image with the retrieval set. Finally, our system warps the annotation results of superpixel matching, and integrates multiple cues in a Markov Random Field (MRF) to obtain an accurate segmentation of the query image. Our non-parametric scene parsing achieves promising results on the LabelMe Outdoor dataset. The system has limited parameters, and captures contextual information naturally in the retrieval and alignment procedure.

### Please cite the following papers if code or part of the code is used :

@inproceedings{naosekpam2019superpixel,
  title={Superpixel Correspondence for Non-parametric Scene Parsing of Natural Images},
  author={Naosekpam, Veronica and Bhowmick, Alexy and Hazarika, Shyamanta M},
  booktitle={International Conference on Pattern Recognition and Machine Intelligence},
  pages={614--622},
  year={2019},
  organization={Springer}
} 

or

Naosekpam, Veronica, Alexy Bhowmick, and Shyamanta M. Hazarika. "Superpixel Correspondence for Non-parametric Scene Parsing of Natural Images." International Conference on Pattern Recognition and Machine Intelligence. Springer, Cham, 2019.

#### And

@inproceedings{naosekpam2019dense,
title={Dense and Partial Correspondence in Non-parametric Scene Parsing},
  author={Naosekpam, Veronica and Paul, Nissi and Bhowmick, Alexy},
  booktitle={International Conference on Machine Intelligence and Signal Processing},
  pages={339--350},
  year={2019},
  organization={Springer}
}

or 

Naosekpam, Veronica, Nissi Paul, and Alexy Bhowmick. "Dense and Partial Correspondence in Non-parametric Scene Parsing." International Conference on Machine Intelligence and Signal Processing. Springer, Singapore, 2019.



