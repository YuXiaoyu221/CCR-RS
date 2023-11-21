# A curvature-driven cloud removal method for remote sensing images

Cloud coverage has become a significant factor affecting the availability of remote sensing images in many applications. To mitigate the adverse impact of cloud coverage and recover ground information obscured by clouds, this paper presents a curvature-driven cloud removal method. Considering that each image can be regarded as a curved surface and the curvature can reflect the texture information well due to its dependence on the surface's undulation degree, the presented method transforms image from natural domain to curvature domain for information reconstruction to maintain details of reference image. In order to improve the overall consistency and continuity of cloud removal results, the optimal boundary for cloud coverage area replacement is determined first to make the boundary pass through pixels with minimum curvature difference. Then the curvature of missing area is reconstructed based on the curvature of reference image, and the reconstructed curvature is inversely transformed to natural domain to obtain a cloud-free image. In addition, considering the possible significant radiometric differences between different images, the initial cloud-free result will be further refined based on specific checkpoints to improve the local accuracy. To evaluate the performance of the proposed method, both simulated experiments and real data experiments are carried out. Experimental results show that the proposed method can achieve satisfactory results in terms of radiometric accuracy and consistency.

if it is used and helpful for your research, please cite:

Xiaoyu Yu, Jun Pan, Mi Wang & Jiangong Xu (2023) A curvature-driven cloud removal method for remote sensing images, Geo-spatial Information Science, DOI: 10.1080/10095020.2023.2189462


Here: "code" gives the matlab code of this algorithm; "data" gives two sets of data for test
