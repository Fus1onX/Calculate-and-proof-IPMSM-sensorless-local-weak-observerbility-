# Observability Analysis of IPMSM

This repository contains MATLAB code for computing the **observability matrix** of Interior Permanent Magnet Synchronous Motors (IPMSMs), its **determinant**, and for verifying the **necessary and sufficient conditions** for rank-based observability criteria.  

The implementation is provided in MATLAB and serves as the supplementary material for the unpublished manuscript **[1]**.  

---

## References

1. *Full-Speed-Range Sensorless Control of IPMSMs With Arbitrary Signal Injection: Observability Analysis and EKF-Based Implementation* (Unpublished manuscript).  
2. P. Vaclavek, P. Blaha, and I. Herman, “AC Drive Observability Analysis,” *IEEE Transactions on Industrial Electronics*, vol. 60, no. 8, pp. 3047–3059, Aug. 2013.  
3. M. Koteich, A. Maloum, G. Duc, and G. Sandou, “Discussion on ‘AC Drive Observability Analysis’,” *IEEE Transactions on Industrial Electronics*, vol. 62, no. 11, pp. 7224–7225, Nov. 2015.  

---

## Repository Structure

The repository contains the following MATLAB scripts:

1. **xxx.m**  
   Reproduces the observability results of IPMSMs discussed in [2] and [3].  
   Includes computation of the observability matrix \( O_{\sigma_1} \) and the determinant terms \( D_{1234} \).  

2. **xxx.m**  
   Provides the proof of **Theorem 1** in Section III-A of [1].  

3. **xxx.m**  
   Provides the proof of **Theorem 2** in Section III-A of [1].  

4. **xxx.m**  
   Provides the proof of **Theorem 3** in Section III-B of [1].  

---

## Usage

1. Clone this repository:
   ```bash
   git clone https://github.com/<your-username>/<repo-name>.git
