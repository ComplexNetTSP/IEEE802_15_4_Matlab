Matlab Code of the implementation of IEEE 802.15.4 Mac/Phy Layer
================================================================

:Info: See `the site <http://bit.ly/T4INBK>`_ for more information. See `GitHub <http://bit.ly/TuTdup>`_ for the latest source.
:Author:Mohamed-Haykel Zayani & Vincent Gauthier
:Maintainer: Vincent Gauthier <vgauthier@luxbulb.org>

Description
-----------
Matlab Code of the implementation of IEEE 802.15.4 Mac/Phy Layer

Introduction
------------
The m-file model an Institute of Electrical and
Electronics Engineers (IEEE) 802.15.4 Machine Access Control (MAC) layer channel
in which multiple non-saturated stations compete to get access to the channel
and to communicate with a sink. It is inspired from the IEEE 802.11 Model [1]
developed by David Griffith and Michael Souryal (Emerging and Mobile Network
Technologies Group, Information Technology Laboratory, National Institute of
Standards and Technology). The objective behind this “adaptation” is to model a
physical layer (PHY), including path loss and shadowing effects which interacts
with the MAC model. The particularity of the proposed model lies in overstepping
the node range disk shaped and taking into consideration the called
“transitional area” [2]_ [3]_. The model relies on the approach of Zuniga and
Krishnamachari [2]_ [3]_. The function ZunPhyModel performs the calculations at the
PHY level and determines the probability of good frame reception towards channel
(signal-to-noise ratio) and radio (modulation and coding) parameters.

References
----------
.. [1] SGIP NIST Smart Grid Collaboration Site, PAP02: Wireless Communicattions for the Smart Grid (6.1.5), Call for Input to Task6, 2010.
.. [2] Zuniga, M. and Krishnamachari, B., “Analyzing the transitional region in low power wireless links,” In Proceedings of the IEEE 1st Annual Conference on Sensor and Ad Hoc Communications and Networks (SECON), pp. 517–526, 2004.
.. [3] Zuniga, M., and Krishnamachari, B., “An analysis of unreliability and asymmetry in lowpower wireless links,” ACM Trans. Sensor Netw., 2007.
