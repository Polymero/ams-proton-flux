/** \class NtpCompact
Compact tree
*/
class NtpCompact {

 public:

  unsigned int status;                      ///< nParticle()+nAntiCluster()*10+nBetaH()*100+nTrTrack()*1000+nTrRecHit()*10000+nTrdCluster()*1000000+nTofClusterH()*100000000
  short int    sublvl1;                     ///< Pattern of LVL1 sub-triggers (8 bits)
  short int    trigpatt;                    ///< Pattern of trigger system members (16 bits)

  float        tof_beta;                    ///< TOF Beta
  float        tof_chisqtn;                 ///< TOF Time Normalized Chi2
  float        tof_chisqcn;                 ///< TOF Spatial Normalized Chi2
  float        tof_q_lay[4];                ///< TOF Charge of each layer
  float        tof_edep_frac[4];            ///< TOF E_{Dep} / Total E_{Dep} on each layer
  short int    tof_trk_ncl;                 ///< Number of TofClusterH matching with associated TrTrack
  short int    tof_beta_patt;               ///< Beta pattern (bit0-3 on: excluded layer, bit4 on: bad beta, bit5 on: bad match TOF/track)
  short int    tof_beta_ncl;                ///< Number of TofClusterH in beta computation
  short int    tof_z;                       ///< Integer charge
  float        tof_z_like;                  ///< -Loglikelihood for integer charge
  short int    tof_z_nhit;                  ///< Number of hits in integer charge computation
  short int    tof_clsn[4];                 ///< [in_time/of_time(2top+2bot)],clusters number

  float        tof_evgeni_beta;             ///< Beta from BetaR

  short int    trk_patty;                   ///< Pattern Y
  short int    trk_pattxy;                  ///< Pattern XY
  float        trk_q_inn;                   ///< Inner Tracker Charge
  float        trk_q_inn_rms;               ///< Inner Tracker Charge RMS
  float        trk_q_lay[9];                ///< Tracker Layer charge (inverted sign for bad status)
  float        trk_edep_frac[9][2];         ///< Tracker on-track E_{dep} over E_{dep} integrated inside 10 cm around the track on each layer/side
  short int    trk_fiducial[5];             ///< Choutko fit pattern in tracker (FS, L1+Inner, L9+Inner, Inner, Inner no MS)
  float        trk_rig[5];                  ///< Choutko fit rigidities (FS, L1+Inner, L9+Inner, Inner, Inner no MS) [GV]
  float        trk_chisqn[5][2];            ///< Choutko fit normalized Chi2 (FS Kalman, FS, L1+Inner, L9+Inner, Inner, Inner no MS | X, Y)
  float        trk_stoermer[5];             ///< Choutko fit directional cutoff for positive particles (FS, L1+Inner, L9+Inner, Inner, Inner no MS) [GV]
  short int    trk_kal_fiducial;            ///< Full-span Kalman fit pattern in tracker 
  float        trk_kal_rig[3];              ///< Full-span Kalman fit rigidity (195, 0, -70 cm) [GV]
  float        trk_kal_chisqn[2];           ///< Full-span Kalman fit normalized Chi2 (X,Y)
  float        trk_kal_stoermer;            ///< Full-span Kalman fit directional cutoff for positive particles [GV]
  float        trk_int_rich_rad[2];         ///< Choutko FS fit interpolation on the RICH radiator (X,Y) [cm] @ Analysis::rich_radiator_z

  float        trk_rho_up;                  ///< R_{Inn-NoMS,up}/R_{Inn-NoMS}
  float        trk_rho_dw;                  ///< R_{Inn-NoMS,down}/R_{Inn-NoMS}

  short int    rich_select;                 ///< Javier selection (see Tools::RichQC)
  unsigned int rich_status;                 ///< nRichHit() + 100*pRichRing->getPMTs() + 10000*pRichRing->getUsedHits()
  float        rich_beta;                   ///< RICH beta best estimator
  float        rich_tot_np;                 ///< Uncorrected total number of p.e. discarding those on PMTs crossed by charged particles
  float        rich_np;                     ///< Uncorrected number of p.e. collected in the ring
  float        rich_np_exp;                 ///< Uncorrected expected number of p.e.
  float        rich_prob;                   ///< Kolmogorov test to the distribution of charge along the ring
  float        rich_bdt;                    ///< RICH BDT from F. Dimiccoli

  float        sa_tof_beta;                 ///< Standalone TOF Beta
  float        sa_tof_chisqtn;              ///< Standalone TOF Time Normalized Chi2
  float        sa_tof_q_lay[4];             ///< Standalone TOF Charge of each layer
  short int    sa_tof_beta_ncl;             ///< Standalone number of TofClusterH in beta reconstruction
  short int    sa_tof_build;                ///< Standalone build number
  short int    sa_tof_clsn[4];              ///< [in_time/of_time(2top+2bot)],clusters number
  float        sa_tof_edep_frac[4];         ///< TOF E_{Dep} / Total E_{Dep} on each layer

  bool         sa_trd_same;                 ///< Is the TRD associated to particle the same of the standalone search? 
  float        sa_trd_q;                    ///< Standalone TRD Charge
  float        sa_trd_chi2;                 ///< Standalone TRD Chi2
  short int    sa_trd_fiducial;             ///< Standalone TRD pattern in tracker
  float        sa_trd_stoermer;             ///< Standalone TRD directional cutoff for positive particles [GV]

  float        sa_exthit_ql1;               ///< Standalone Unbiased L1 Hit Charge
  float        sa_exthit_dl1[2];            ///< Standalone Unbiased L1 Hit Distance from TRD Track (X-closest, Y)
  short int    sa_exthit_status_l1;         ///< Standalone Unbiased L1 Hit Status (L1)

  float        sa_ecal_edepd;               ///< Standalone Largest ECAL Shower EDep [GeV]

  float        mc_momentum;                 ///< MC generated momentum 
