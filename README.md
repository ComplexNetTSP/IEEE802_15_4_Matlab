IEEE802_15_4_Matlab
===================
Mohamed-Haykel Zayani & Vincent Gauthier 
----------------------------------------

Laboratory CNRS S.A.M.O.V.A.R. - Dept RS2M
------------------------------------------

Description 
-----------
Matlab Code of the implementation of IEEE 802.15.4 Mac/Phy Layer

Introduction
------------
The m-file model an Institute of Electrical and Electronics Engineers (IEEE) 802.15.4 Machine Access Control (MAC) layer channel in which multiple non-saturated stations compete to get access to the channel and to communicate with a sink. It is inspired from the IEEE 802.11 Model [1] developed by David Griffith and Michael Souryal (Emerging and Mobile Network Technologies Group, Information Technology Laboratory, National Institute of Standards and Technology). The objective behind this “adaptation” is to model a physical layer (PHY), including path loss and shadowing effects which interacts with the MAC model. The particularity of the proposed model lies in overstepping the node range disk shaped and taking into consideration the called “transitional area” [2, 3]. The model relies on the approach of Zuniga and Krishnamachari [2, 3]. The function ZunPhyModel performs the calculations at the PHY level and determines the probability of good frame reception towards channel (signal-to-noise ratio) and radio (modulation and coding) parameters.
