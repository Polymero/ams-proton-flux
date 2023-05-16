#ifndef __Ntp_h__
#define __Ntp_h__

#include "TObject.h"
#include "TString.h"

#include <cmath>
#include <algorithm>
#include <string>

using namespace std;

class RichBDTMgr;
class TrackerBDTMgr;
class ClassifierManager;
class RichOccupancy;

/** \class RTIInfo
Real Time Information, calculated and stored for each second.
It is a full copy of the AMSSetupR::RTI database.
*/
class RTIInfo {

 public:

  bool         isinsaa;     ///< Is Inside SAA
  unsigned int run;         ///< Run
  int          evno;        ///< First event number in the RTI second
  int          evnol;       ///< Last event number in the RTI second
  float        lf;          ///< Livetime [0,1]
  int          ret_cf[7];   ///< Max geomegnetic cutoff calculation status (>=0: success, -1: failure, <=-2: error)
  float        cf[7][4][2]; ///< Max geomagnetic cutoff in the field of view (Stoermer,IGRF-RTI,IGRF-12,IGRF-12+,IGRF-12++,Tsy05(Max-Sec),Tsy05(Min-Pri)|25,30,35,40 degrees|-,+) [GV]
  float        mphe;        ///< most probable He rigidity
  float        theta;       ///< Theta GTOD [rad]
  float        phi;         ///< Phi GTOD [rad]
  float        r;           ///< Altitude GTOD [cm]
  float        zenith;      ///< AMS inclination w/ the zenith [degrees]
  float        glat;        ///< Pointing galactic latitude, -1 if failed [degrees]
  float        glong;       ///< Pointing galactic longitude, -1 if failed [degrees]
  float        nev;         ///< Number of events on disk (nev+nerr = all events)
  float        nerr;        ///< Number of missing events
  float        ntrig;       ///< Number of events with trigger
  float        nhwerr;      ///< Number of events with DAQ error (from JINJstatus)
  float        npart;       ///< Number of events with a particle with tof+tracker+rich+ecal
  float        nl1l9[2][2]; ///< Events with L1 and/or L9 hits (L1,L9|X,Y)
  float        dl1l9[2][3]; ///< Mean difference bewteen PG ad CIEMAT alignment of L1 and L9 (L1,L9|X,Y,Z) [um]
  float        mtrdh;       ///< Average number of TrdRawHit per event

  //! 0 if good, otherwise bad
  /** bitcode:
    - bit0: duplicated events
    - bit1: event number flip
    - bit2: event missing at the beginging of second
    - bit3: event missing at the end of second
    - bit4: second at the begining of run
    - bit5: second at the end of run
   */
  int          good;
  unsigned int utime;       ///< JMDC unix time [s]
  unsigned int usec[2];     ///< JMDC unix time microsecond for first and last event [us]
  double       utctime[2];  ///< UTC time for first and last event [s]
  double       betasun;     ///< solar beta angle [degree]

  //! The standard RTI selection for fluxes determination (from Q. Yan)
  bool Select();
  //! A list of bad runs (from S. Haino)
  bool IsBadRun();

  RTIInfo(){}
  virtual ~RTIInfo(){}
  ClassDef(RTIInfo,1);
};

/** \class FileInfo
Informations related to a single processed AMS ROOT file (or to a part of it).
*/
class FileInfo {

 public:

  unsigned int run;         ///< Run
  unsigned int utime[2];    ///< JMDC unix time of first and last event [s]
  int          event[2];    ///< Event number of first and last event
  int          nentries;    ///< Number of entries

  FileInfo(){}
  virtual ~FileInfo(){}
  ClassDef(FileInfo,1);
};

/** \class FileMCInfo
Informations related to a single processed AMS ROOT MC file (only stored in case of MC processing).
Includes MC generator informations (min/max momentum, number of generated events, ...).
*/
class FileMCInfo {

 public:

  int   event[2];           ///< Min and max event number (event[1]-event[0]~ngen)
  int   charge;             ///< Primary charge [e]
  float mass;               ///< Primary mass [GeV/c^2]
  float momentum[2];        ///< Min and max momentum [GeV/c]
  int   pid_datacard;       ///< Particle ID from datacard (G3 PID, if not available PDG PID)
  float mom_datacard[2];    ///< Min and max momentum from datacard [GeV/c]
  int   ngen_datacard;      ///< Number of generated events from datacard
  float mom_filename[2];    ///< Min and max momentum from file [GeV/c]
  bool  isl1_filename;      ///< L1 focus
  bool  isl19_filename;     ///< L1&L9 focus
  bool  istb_filename;      ///< Test beam MC

  //! best guess of number of generated
  int GetNGen();
  //! best guess of minimum rigidity generated
  double GetRMin();
  //! best guess of maximum rigidity generated
  double GetRMax();

  FileMCInfo(){}
  virtual ~FileMCInfo(){}
  ClassDef(FileMCInfo,1);
};

/** \class ProcInfo
Information realted to the data reduction program execution (execution time, executable version, ...).
There is one and only one of these for each ntuple produced.
*/
class ProcInfo {

 public:

  unsigned int utctime[2];    ///< Processing UBegin/UTerminate UTC time
  int          nevents[6];    ///< Events selection (processed,good daq/trig,have beta and/or track,after prescaling,good refit,filling)
  float        time[15];      ///< Processing times (Init,Reading,Select/Prescale,Clear,Header,HeaderMC,TRD,TOF,Tracker,RICH,ECAL,Anti,Stand-alone,Fill)
  bool         is_mc;         ///< Is MC
  int          pres_strategy; ///< Prescaling strategy used
  int          lib_ver[2];    ///< Version and subversion of the ntuple library
  int          exe_ver[2];    ///< Version and subversion of the executable
  char         ams_ver[100];  ///< AMS software version

  ProcInfo(){}
  virtual ~ProcInfo(){}
  ClassDef(ProcInfo,1);
};

/** \class NtpSHeader
Short header information with run, event and "RTI time". Separated to have fast index building.
*/
class NtpSHeader {

 public:

  unsigned int run;    //< Run
  int          event;  //< Event
  unsigned int utime;  //< JMDC unix time [s]
  unsigned int herror; //< AMS Header error

  NtpSHeader(){}
  virtual ~NtpSHeader(){}
  ClassDef(NtpSHeader,1);
};

/** \class NtpHeader
Informations about prescaling, reconstruction, data acquisition, trigger, high-level infos and orbit.
*/
class NtpHeader {

 public:

  float    pres_weight;       ///< Prescaling weight
  float    pres_trck_rig;     ///< Tracker rigidity used in precaling
  float    pres_tofh_beta;    ///< TOF beta used in prescaling
  float    pres_rich_beta;    ///< RICH beta used in prescaling
  float    pres_trck_inn_q;   ///< Tracker charge used in prescaling
  long int pres_patt;         ///< Pattern of conditionals used for prescaling
  double   utc_time;          ///< UTC time [s]
  double   utc_time_error;    ///< UTC error [s]

  int      nparticle;         ///< Number of Particle
  int      nmceventg;         ///< Number of MCEventg
  int      nanti;             ///< Number of AntiClusteir
  int      ntof;              ///< Number of TofCluster
  int      ntofh;             ///< Number of TofClusterH
  int      ntrtrack;          ///< Number of TrTrack
  int      ntrrechit;         ///< Number of TrRecHit
  int      ntrdtrack;         ///< Number of TrdTrack
  int      ntrdhtrack;        ///< Number of TrdHTrack
  int      ntrdsegment[2];    ///< Number of TrdSegment (X,Y)
  int      ntrdcluster;       ///< Number of TrdCluster
  int      nrich;             ///< Number of RichRing
  int      nrichhit;          ///< Number of RichHit
  int      nrichhit_crossed;  ///< Number of crossed RichHit
  int      necal;             ///< Number of EcalShower

  //! Event error
  /** bitcode:
    - bit0: error on bits (0x7F00)c on the four JINJStatus word
    - bit1: whatever bad REPLY at JINJ, JINF or LVL1 (ROOM, Sync, 0-replies, ...)
    - bit2: desync error at whatever level
    - bit3: decoding procedure error
    - bit4: whatever board with bad REPLY
    - bit5: reconstruction error
   */
  int      error;
  int      size;              ///< Event length in bytes
  float    livetime;          ///< Livetime
  int      antipatt;          ///< ACC sectors in coincidence with FT (8 bits)
  int      sublvl1;           ///< Pattern of LVL1 sub-triggers (8 bits)
  int      trigpatt;          ///< Pattern of trigger system members (16 bits)
  float    rads;              ///< ISS orbit altitude (GTOD, prediction from NORAD web site, accuracy ~1 km) [cm]
  float    thetas;            ///< ISS theta (GTOD, prediction from NORAD web site, accuracy ~1 km) [rad]
  float    phis;              ///< ISS phi (GTOD, prediction from NORAD web site, accuracy ~1 km) [rad]
  float    yaw;               ///< ISS yaw (LVLH) [rad]
  float    pitch;             ///< ISS picth (LVLH) [rad]
  float    roll;              ///< ISS roll (LVLH) [rad]
  float    velocitys;         ///< ISS velocity [rad/s]
  float    velthetas;         ///< ISS velocity theta [rad]
  float    velphis;           ///< ISS velocity phi [rad]
  float    thetam;            ///< Magnetic theta (eccentric dipole model) [rad]
  float    phim;              ///< Magnetic phi (eccentric dipole model) [rad]

  float    dedx_beta[4];      ///< Beta estimation from dE/dx (Tracker,TOF,TRD,Tracker+TOF+TRD)
  float    mass_quality;      ///< Mass quality -Log(Likelihood) estimator
  int      ntrdseg_vertex;    ///< Number of TRD segments which make a vertex above
  float    bl2;               ///< Int(BL^2) on the trajectory between L2 and L7 with step of 1 cm [T m^2]


  //! The unbiased-charged trigger subdetector flag (TOF 3/4)
  bool IsFTCP0() { return (trigpatt&0x2)!=0; }
  //! Charged particle trigger for physics analysis
  bool IsChargedPhysTrigger() { return IsFTCP0()&&((sublvl1&0x3E)!=0); }
  //! Charged particle trigger for trigger efficiency calculation (to be added with factor 100 to trigger efficiency denominator)
  bool IsChargedUnphysTrigger() { return IsFTCP0()&&((sublvl1&0x3E)==0); }

  NtpHeader(){}
  virtual ~NtpHeader(){}
  ClassDef(NtpHeader,2);
};

/** \class NtpMCHeader
Informations related to the MC truth (only stored in case of MC processing).
*/
class NtpMCHeader {

 public:

  int   pid;              ///< Particle ID (GEANT3 ID if existing, else PDG ID)
  float mass;             ///< Mass [GeV/c^2]
  float charge;           ///< Charge
  float momentum[22];     ///< Momentum (at 22 different heights, if surviving)
  float coo[22][3];       ///< Coordinates (at 22 different heights, if surviving) [cm]
  float dir[3];           ///< Generator direction

  int   hit_pid[9];       ///< Particle ID (GEANT3 ID if existing, else PDG ID) of a matched TrMCCluster along the first track (=0 if no cluster is found)
  float hit_pos[9][3];    ///< Particle postion a tracker layer
  float hit_mom[9][3];    ///< Particle postion a tracker layer
  float hit_edep[9];      ///< Particle energy deposition

  int   rich_part[2][10]; ///< PID of first 10 particles generating detected photons on the RICH (in radiator, in the PMT plane)
  float rich_mom[2][10];  ///< Momentum of first 10 particles generating detected photons on the RICH (in radiator, in the PMT plane) [GeV]
  int   rich_np[2][10];   ///< Number of photons detected photons on the RICH by the first 10 particles (in radiator, in the PMT plane)

  NtpMCHeader(){}
  virtual ~NtpMCHeader(){}
  ClassDef(NtpMCHeader,2);
};

/** \class NtpTrd
Informations related to the TRD reconstruction.
Variables from three different approaches are stored: 1) from the orginal TrdTrack reconstruction (MIT),
and 2) from the TrdK reconstruction (MIT).
*/
class NtpTrd {

 public:

  float theta;                 ///< Theta [rad]
  float phi;                   ///< Phi [rad]
  float coo[3];                ///< Track passing point [cm]
  float chisq;                 ///< Chi2/NDF of a linear fit
  int   nseg;                  ///< Number of segments in the TrdTrack
  float q;                     ///< Charge
  bool  trdk_is_align_ok;      ///< TrdK alignment is loaded
  bool  trdk_is_calib_ok;      ///< TrdK calibration is loaded
  float trdk_q[3];             ///< TrdK charge (dE/dx&&delta-rays,dE/dx,delta-rays)
  float trdk_q_err[3];         ///< TrdK charge error
  float trdk_q_upper;          ///< TrdK charge upper
  float trdk_q_lower;          ///< TrdK charge lower
  int   trdk_q_nhit;           ///< TrdK number of hits used for charge estimation
  int   trdk_q_nhit_refit;     ///< TrdK number of hits used for charge estimation refit
  int   trdk_q_nhit_nuclei;    ///< TrdK number of hits used for charge estimation dE/dx
  int   trdk_q_nhit_dr;        ///< TrdK number of hits used for charge estimation delta-rays
  int   trdk_nhit;             ///< TrdK number of hits used
  float trdk_ampl[25];         ///< Energy deposit, pressure-corrected [ADC] (up to 40 hits would be possible)
  float trdk_path[25];         ///< Pathlengths in cm (up to 40 hits would be possible)
  int   trdk_like_valid[3];    ///< TrdK valid likelihood calculation (inner,max span,refit)
  float trdk_like_e[3];        ///< TrdK -loglikelihood for electron (inner,max span,refit)
  float trdk_like_p[3];        ///< TrdK -loglikelihood for proton (inner,max span,refit)
  float trdk_like_he[3];       ///< TrdK -loglikelihood for helium (inner,max span,refit)
  float trdk_like_d[2];        ///< TrdK -loglikelihood for deuteron (inner,max span)
  int   trdk_like_nhit[3];     ///< TrdK number of hits used for likelihood calculation (inner,max span,refit)
  int   trdk_like_nhit_off[3]; ///< TrdK number of hits off the likelihood calculation (inner,max span,refit)
  float trdk_like_ampl_off[3]; ///< TrdK total amplitude of hits off the likelihood calculation (inner,max span,refit)
  int   vertex_nx;             ///< Number of reconstructed vertexes along XZ projection
  int   vertex_ny;             ///< Number of reconstructed vertexes along YZ projection
  int   vertex_nsegx;          ///< Number of segments in best 3D rec. vertex along XZ projection
  int   vertex_nsegy;          ///< Number of segments in best 3D rec. vertex along YZ projection
  float vertex_coo[3];         ///< Coordinate of best 3D rec. vertex (x,y,z) [cm]
  float vertex_d2;             ///< Summed squared distance of 3D rec. vertex from XZ and YZ segments (used for minimization) [cm^2]

  //! linear interpolation
  void Interpolate(float z, float xy[2]);
  //! is inside RICH radiator (0: no, 1: inside NaF, 2: inside Aerogel)
  int IsInsideRich();
  //! Tracker fiducial volume pattern (bit0: is in L1, bit1: is in L2, ... bit8: is in L9)
  int GetPatternInsideTracker();
  //! is inside TRD
  bool IsInsideTRD();
  //! is inside ECAL
  bool IsInsideECAL();

  //! likelihood ratio as defined in TrdK (-log(p/(p+He)))
  double trdk_like_phe(int i);
  //! likelihood ratio as defined in TrdK (-log(e/(e+p)))
  double trdk_like_ep(int i);

  NtpTrd(){}
  virtual ~NtpTrd(){}
  ClassDef(NtpTrd,1);
};

/** \class NtpTof
Informations related to the TOF reconstruction (from both Q.Yan and E.Choumilov reconstructions).
*/
class NtpTof {

 public:

  float beta;              ///< Beta
  float beta_err;          ///< Beta fitting error
  float t0;                ///< Time offset
  int   beta_ncl;          ///< Number of TofClusterH in beta computation
  int   beta_patt;         ///< Beta pattern (0: 4 measurement, !=0: each bit indicates which layer is excluded)
  float chisqcn;           ///< Normalized spatial Chi2/NDF with respect to associated TrTrack
  float chisqtn;           ///< Normalized time fitting Chi2/NDF
  float theta;             ///< TOF Track theta [rad]
  float phi;               ///< TOF Track phi [rad]
  float coo[3];            ///< TOF Track interpolation at z=0 [cm]
  float q;                 ///< Charge
  float q_err;             ///< Charge error
  int   q_nhit;            ///< Number of hits in charge computation
  int   z;                 ///< Integer charge
  float z_like;            ///< -Loglikelihood for integer charge
  int   z_nhit;            ///< Number of hits in integer charge computation
  float time[4];           ///< Layer time
  float q_lay[4];          ///< Layer charge
  float q_lay_uncorr[4];   ///< Layer charge with no beta/rigidity correction
  int   flag;              ///< Beta flag, bit0:!IsGoodBeta, bit1:!IsTkTofMatch
  int   trk_ncl;           ///< Number of TofClusterH matching with associated TrTrack

  //! Plane flag (<0 no cluster)
  /** The flag is defined as A+100*B+10000*C, where:
    - A is for the counter:
      -# bit0: !IsIsolation (no hit in the neightboring counters)
      -# bit1: IsOnOverlap (the particle is passing within 1 cm from the side of the counter)
    - B is for side 0 electronics:
      -# bit0: !TestExistHS (1 if there is no TofRawSide)
      -# bit1: !IsGoodSide (1 if no FT signal back from LVL1 | no LT (time meas.) signal in a good time window from FT | no ADC from anode)
      -# bit2: !IsOneLT (1 if more than one LT in a good time window from FT)
      -# bit3: !IsExistHT (1 if no HT trigger signal, in the cumulative HT channel)
    - C is for side 1 electronics:
      -# bit0: !TestExistHS (1 if there is no TofRawSide)
      -# bit1: !IsGoodSide (1 if no FT signal back from LVL1 | no LT (time meas.) signal in a good time window from FT | no ADC from anode)
      -# bit2: !IsOneLT
      -# bit3: !IsExistHT
  */
  int   flagp[4];
  int   nclp[4][3];        ///< Number of TofClusterH (plane=0,1,2,3|0:used,1:nearby,2:others)
  float edep[4][5];        ///< Energy deposited (plane=0,1,2,3|0:used,1:nearby,2:others,3:max not nearby,4:overlap not in BetaH)
  float ovresm[4];         ///< resmeasoverlap (?)
  float ovresb[4];         ///< resborderoverlap (?)
  float clsdt[4];          ///< [in_time/of_time(2top+2bot)],aver.times wrt beta-used hit
  float clsed[4];          ///< [in_time/of_time(2top+2bot)],tot.edep
  short int clsn[4];       ///< [in_time/of_time(2top+2bot)],clusters number
  float evgeni_beta;       ///< Beta from BetaR
  int   evgeni_beta_patt;  ///< Beta pattern from BetaR

  //! Calculate the average charge combining several TOF layers
  double GetQ(int &n, double& rms, int pattern, bool good_flag = false, bool uncorr = false);

  NtpTof(){}
  virtual ~NtpTof(){}
  ClassDef(NtpTof,1);
};

/** \class NtpTracker
Informations realted to the Tracker reconstruction.

Most of the informations are relative to the TrTrack selected in Analysis.cxx.
In particular we included several refits with different algorithms and different patterns.

Information about the highest momentum track not used in the analysis are
also stored, as well as informations about cluster morphology and isolation.
*/
class NtpTracker {

 public:

  int    patty;                 ///< Bit pattern of Y clusters
  int    pattxy;                ///< Bit pattern of XY hits
  int    s0;                    ///< Z=0 sensor grid sequence
  int    s2;                    ///< L2 sensor grid sequence
  int    s8;                    ///< L8 sensor grid sequence
  float  q_inn[3];              ///< Inner tracker charge (Old, Y. Jia, H. Liu)
  int    q_inn_nhit[3];         ///< Number of clusters used for inner tracker charge calculation (Old, Y. Jia, H. Liu)
  float  q_inn_rms[3];          ///< Inner tracker charge error (Old, Y. Jia, H. Liu)
  float  q_lay[3][9];           ///< Layer charge X/Y-side (Old, Y. Jia, H. Liu)
  float  q_lay_uncorr[3][9];    ///< Layer charge X/Y-side (Old, Y. Jia, H. Liu, no beta corr.)
  float  q_clu[2][9];           ///< Cluster charge X/Y-side (old algorithm)
  float  q_clu_uncorr[2][9];    ///< Cluster charge X/Y-side (old algorithm, no beta corr.)
  int    q_clu_status[2][9];    ///< Cluster charge X/Y-side status (old algorithm)
  //! Refit rigidity (interpolation @ z=195,0 and -70 cm) [GV]
  /* - 0: Max available span, Kalman fit (algo=6)
     - 1: Inner Tracker, Kalman fit (algo=6)
     - 2: Upper Inner Tracker (exclude Layer 7 and 8), Kalman fit (algo=6)
     - 3: Lower Inner Tracker (exclude Layer 2), Kalman fit (algo=6)
     - 4: Inner Tracker + Layer 1, Kalman fit (algo=6)
     - 5: Inner Tracker + Layer 9, Kalman fit (algo=6)
     - 6: Inner Tracker, Vitaly fit (algo=1)
     - 7: Inner Tracker, Vitaly fit no MS (algo=21)
     - 8: Inner Tracker, Chickanian fit (algo=3)
  */
  float  rig[9][3];
  float  invrigerr[9];          ///< Refit inverse rigidity error [1/GV]
  float  theta[9];              ///< Refit theta [rad]
  float  phi[9];                ///< Refit phi [rad]
  float  coo[9][3];             ///< Refit track passing point [cm]
  float  chisqn[9][2];          ///< Refit Chi2/NDF of X/Y-side
  double int_l2[9][2];          ///< Refit L2 interpolation (X,Y) [cm] @ Analysis::tracker_layers_z[1]
  double coo_lay[9][3];         ///< TrRecHit coordinates [cm]
  double int_lay[9][3];         ///< Reference fit interpolation on each layer (X,Y,Z) [cm] @ Analysis::tracker_layers_z[]
  float  int_rich_rad[2];       ///< Reference fit interpolation on the RICH radiator (X,Y) [cm] @ Analysis::rich_radiator_z
  float  theta_rich_rad;        ///< Reference fit interpolation on the RICH radiator theta [rad] @ Analysis::rich_radiator_z
  float  phi_rich_rad;          ///< Reference fit interpolation on the RICH radiator phi [rad] Analysis::rich_radiator_z
  float  int_rich_pmt[2];       ///< Reference fit interpolation on the RICH PMT plane (X,Y) [cm] @ Analysis::rich_pmt_plane_z
  float  stoermer_cutoff[2];    ///< Reference fit directional Stoermer cutoff (+,-) [GV]
  float  igrf_cutoff[2];        ///< Reference fit directional IGRF cutoff (+,-) [GV]
  float  igrf_cutoff_upper[2];  ///< Reference fit directional IGRF Upper cutoff (+,-) [GV]
  float  igrf_cutoff_lower[2];  ///< Reference fit directional IGRF Lower cutoff (+,-) [GV]
  float  inn_rig[7];            ///< Inner tracker refit rigidity (refit excluding layers one-by-one) [GV]
  float  inn_invrigerr[7];      ///< Inner tracker refit inverse rigidity error (refit excluding layers one-by-one) [1/GV]
  float  inn_res[7][2];         ///< Inner tracker refit residuals for X/Y-side (refit excluding layers one-by-one) [cm]
  float  inn_chisqn[7][2];      ///< Inner tracker refit Chi2/NDF of X/Y-side (refit excluding layers one-by-one)
  int    nclu_lay[2][9][7];     ///< Total number of clusters (X,Y|layer|on-track,0.1cm,1cm,2cm,5cm,10cm,all layer)
  float  edep_lay[2][9][7];     ///< Total energy deposition (X,Y|layer|on-track,0.1cm,1cm,2cm,5cm,10cm,all layer) [keV]
  float  max_edep_lay[2][9];    ///< Maximum energy deposition on a layer (X,Y|layer) [keV]
  float  d2_max_edep_lay[2][9]; ///< Squared distance of the maximum energy deposition cluster from the track (X,Y|layer) [cm^2]
  float  signal_ratio[7];       ///< Ratio between cluster amplitude and neightboring 10 strips
  float  feet_dist[7];          ///< Distance from the closet ladder feet
  int    sec_patty;             ///< Second track bit pattern of Y cluster
  int    sec_pattxy;            ///< Second track bit pattern of XY hits
  float  sec_inn_q;             ///< Second track inner tracker charge (no rig. corr.)
  float  sec_inn_rig;           ///< Second track inner tracker rigidity [GV]
  float  sec_inn_chisqn[2];     ///< Second track inner tracker Chi2/NDF of X/Y-side

  float  scat_rigi_theta[2][2];       ///< Scattering angle using only track informations (materials between L34,L56|x,y) [degrees]
  float  scat_beta_theta[2][2][2][2]; ///< Scattering angle using track and velocity (beta TOF,RICH|mass p,d|materials between L34,L56|x,y) [degrees]

  //! is inside RICH radiator (0: no, 1: inside NaF, 2: inside Aerogel)
  int IsInsideRich();
  //! Tracker fiducial volume pattern (bit0: is in L1, bit1: is in L2, ... bit8: is in L9)
  int GetPatternInsideTracker();
  //! is inside TRD
  bool IsInsideTRD();
  //! is inside ECAL
  bool IsInsideECAL();

  NtpTracker(){}
  virtual ~NtpTracker(){}
  ClassDef(NtpTracker,3);
};

/** \class NtpRich
Informations related the RICH reconstruction (both CIEMAT and LIP reconstructions are included).

In case of no ring reconstructed an estimation of the number of photo-electrons in the case of electrons
or protons (with rigidity from Tracker) is given, with the purpose of creating a veto for fast particles.
*/
class NtpRich {

 public:

  // global

  int   nparticle;                ///< Number of crossed PMTs (= number of traversing particles)
  int   tot_hit_uncorr;           ///< Uncorrected total number of hits
  int   tot_pmt_uncorr;           ///< Uncorrected total number of PMTs
  float tot_p_uncorr;             ///< Uncorrected total number of p.e. discarding those on PMTs crossed by charged particles
  float max_p_uncorr[2];          ///< Uncorrected maximum number of p.e. in a PMT (all, excluding crossed PMTs)

  int   pmt_nhit_uncorr[5];       ///< First 5 PMTs by number of p.e., uncorrected number of hits
  float pmt_np_uncorr[5];         ///< First 5 PMTs by number of p.e., uncorrected number of p.e.
  float pmt_dist[5];              ///< First 5 PMTs by number of p.e., CoG distance from reference track [cm]
  int   pmt_pmt[5];               ///< First 5 PMTs by number of p.e., PMT number

  int   tot_hit[2][5];            ///< Number of hits (0:In-Ring, 1:Out-of-Ring | 0:Bad, 1:Crossed by Primary, 2:Crossed by Secondary, 3:Carlos-Crossed, 4:Good)
  float tot_p[2][5];              ///< Number of photons (0:In-Ring, 1:Out-of-Ring | 0:Bad, 1:Crossed by Primary, 2:Crossed by Secondary, 3:Carlos-Crossed, 4:Good)

  // hits

  //! First 30 hit status (first hits with beta-hyp.)
  /** Status bits:
   - 1: Hit used in the ring number 1
   - 2: Hit used in the ring number 2
   - 3: Hit used in the ring number 3
   - ...
   - 10: Hit used in the ring number 10 (no more than 10 created)
   - 29: Channel taggeg as good in calibration
   - 30: Gain mode chosen for the hit 0=x1(low) 1=x5(high)
   - 31: Hit belongs to a PMT apparently crossed by a charged particle
  */
  unsigned int hit_stat[30];
  unsigned int hit_stat2[30];            ///< First 30 hit additional status (bad, bad occ., primary cross., secondary cross., not in sel. ring)
  int          hit_used[30];             ///< First 30 hit used-word (0: used in ring as direct, 1: used in ring as reflected, otherwise not used in any ring)
  int          hit_chan[30];             ///< First 30 hit channel number (16*PMT+pixel)
  float        hit_np_uncorr[30];        ///< First 30 hit uncorrected number of p.e.
  float        hit_beta[30][2];          ///< First 30 hit raw beta (0: direct, 1: reflected, -1 for bad evaluation)
  int          tot_hyp_hit_uncorr[2][3]; ///< Total number of hit compatible with beta=1 hypothesys (all, out of ring|direct, reflected, reflected)
  float        tot_hyp_p_uncorr[2];      ///< Total number of p.e. for beta=1 hypothesys (all, out of ring)

  // ring

  int   selection;                ///< Javier selection (see Tools::RichQC)
  int   is_naf;                   ///< Used radiatior is NaF
  int   status;                   ///< Status word
  int   correction_status;        ///< Charge corrections fail flag (-1/0/1 : Not/Done/Failed)
  int   occupancy_status;         ///< Number of channels tagged in occupancy check (-1 Failed)
  int   nhit;                     ///< Number of hits in the ring
  int   nhit_uncorr;              ///< Uncorrected number of hits in the ring
  int   nhit_refl;                ///< Number of hits which are consistent with reflected photons
  int   npmt;                     ///< Number of PMTs in the ring
  int   npmt_uncorr;              ///< Uncorrected number of PMTs in the ring
  float np;                       ///< Number of p.e. collected in the ring
  float np_uncorr;                ///< Uncorrected number of p.e. collected in the ring
  float np_w[10];                 ///< Photoelectrons associated to the ring for different windows sizes (1,2,3...10)
  float npmt_exp;                 ///< Expected number of PMTs
  float np_exp;                   ///< Expected number of p.e. for a Z=1 ring with reconstruction input pars of the current event
  float np_exp_uncorr;            ///< Uncorrected expected number of p.e.
  float np_exp_elec;              ///< Expected number of p.e. for a Z=1 ring with reconstruction input pars of the current event and beta = 1
  float prob;                     ///< Kolmogorov test to the distribution of charge along the ring
  float width;                    ///< Width of the distribution of charge around the ring over the expected one
  float udist;                    ///< (1/^2) for unused hits which do not belong to PMTs crossed by a charged particle (?)
  int   tile_id;                  ///< Tile id for the tile crossed by the particle
  float rad_coo[2];               ///< Track position interpolated to the radiator entrance (as used in the reconstruction) [cm]
  float rad_theta;                ///< Track theta interpolated to the radiator entrance (as used in the reconstruction) [rad]
  float rad_phi;                  ///< Track phi interpolated to the radiator entrance (as used in the reconstruction) [rad]
  float distance_tile_border;     ///< Distance of the track impact point in the radiator to the border of the radiator tile
  float beta;                     ///< Beta (tileCorrection)
  float beta_res;                 ///< Expected resolution
  float beta_rms;                 ///< Expected resolution RMS
  float beta_raw;                 ///< Raw beta
  float beta_refit;               ///< Refit beta
  float beta_corrected;           ///< Beta corrected for impact point and direction (best beta estimator)
  float beta_unifcorr;            ///< Beta corrected for impact point and direction (with previous correction, only for comparison)
  float q;                        ///< Charge
  float q_consistency;            ///< Statistical test to check if the hit by hit charge is consistent PMT-by-PMT
  float q_res;                    ///< Expected charge resolution
  float q_rms;                    ///< Expected charge resolution RMS
  int   nclus;                    ///< Number of beta clusters
  int   clus_size[10];            ///< Size of first 10 clusters (ordered by size)
  float clus_mean[10];            ///< Average beta of first 10 clusters (ordered by size)
  float clus_rms[10];             ///< RMS of first 10 clusters (ordered by size)
  float lip_beta;                 ///< LIP beta
  float lip_q;                    ///< LIP charge

  // veto

  float veto_np_exp_elec;         ///< Estimate number of expected p.e. in the case of electron (beta = 1) and no ring present
  float veto_np_exp_prot;         ///< Estimate number of expected p.e. in the case of proton (beta from rigidity) and no ring present
  float veto_beta_trk_prot;       ///< Beta estimated from rigidity

  //! Get total number of photons for a give beta hypothesys
  /** parameters:
    - used: XY,
      - X: 0, all; 1: not used in the ring
      - Y: 0, all; 1: not in crossed PMT (default); 2: not in crossed PMT (mine); 3: not in crossed PMT (both)
    - res[2]: resolution  applied for NaF (0), and Aerogel (1)
  */
  double GetNpBetaHyp(int& nhit, double beta, int used, double* res = 0);

  NtpRich(){}
  virtual ~NtpRich(){}
  ClassDef(NtpRich,1);
};

/** \class NtpEcal
Informations related to the two highest energetic ECAL showers.
The last index of every variable in this class refers to the two shower (the most energetic one first).
No geometrical matching with TrTrack, TrdTrack or BetaH is requested.
*/
class NtpEcal {

 public:

  float energyE[2];          ///< Reconstructed energy (GeV)
  float energyD[2];          ///< Energy deposit (GeV)
  float entry[3][2];         ///< Entry point (X,Y,Z) [cm]
  float exit[3][2];          ///< Exit point (X,Y,Z) [cm]
  float cog[3][2];           ///< Center-of-gravity (X,Y,Z) [cm]
  float dir[3][2];           ///< Direction
  float chi2dir[2];          ///< Chi2 direction
  float mips_tag[2];         ///< Nacho's MIP tag
  float moliere[2];          ///< energy(+-3 cm)/energy ratio
  float rear_leak[2];        ///< Rear leak
  float tmax[2];             ///< Shower maximum (cm)
  float long_disp[2];        ///< Longitudinal dispersion
  int   nhit[2];             ///< Total number of hits with good status
  float bdt[2];              ///< Boosted decision tree classifier
  float q[2];                ///< Charge
  float edep_tot_lay[18][2]; ///< Total energy deposited on a single layer
  float edep_max_lay[18][2]; ///< Maximum energy deposited on a cell of a single layer
  int   nhit_lay[18][2];     ///< Number of hits per layer

  NtpEcal(){}
  virtual ~NtpEcal(){}
  ClassDef(NtpEcal,1);
};

/** \class NtpAnti
Informations realted to anti-counters (PG reconstruction alone).
*/
class NtpAnti {

 public:

  int   npair[8];  ///< number of pairs (if <0 number of time history of the found side)
  float chisq[8];  ///< minimum chisquare respect to the T&Z guess
  float unfz[8];   ///< unfolded zeta in the [-40,40] range less than (6cm resolution 10cm RMS)
  float edep[8];   ///< evaluation of energy deposition (MeV) from rawq

  NtpAnti(){}
  virtual ~NtpAnti(){}
  ClassDef(NtpAnti,1);
};

/** \class NtpStandAlone
Contains informations of a TOF(TRD) track built without relying on the Tracker track.

A sample of TOF tracks independent from Tracker track is obtained running the
TOF reconstruction to obtain BetaH object from TrdTrack, EcalShower or by combination
of TofClusterH.

From the BetaH we derive the association with the TrdTrack, if existing.
Then the TrdKCluster object is created from the TrdTrack, if existing,
otherwise from the BetaH.

Eventually tracker hits on the external layers (L1,L9) with the largest charge
and closest to the TOF(TRD) track are also stored.
*/
class NtpStandAlone {

 public:

  float theta;                    ///< Theta of track (TrdTrack or, if not found, BetaH) [rad]
  float phi;                      ///< Phi of track (TrdTrack or, if not found, BetaH) [rad]
  float coo[3];                   ///< Coordinate of track of track (TrdTrack or, if not found, BetaH) [cm]
  float beta;                     ///< Beta
  float beta_err;                 ///< Beta fitting error
  int   beta_ncl;                 ///< Number of TofClusterH in beta computation
  int   beta_patt;                ///< Beta pattern (0: 4 measurement, !=0: each bit indicates which layer is excluded)
  float beta_chisqt;              ///< Time fitting Chi2 {TO BE REMOVED}
  float beta_chisqtn;             ///< Normalized time fitting Chi2/NDF
  float beta_q;                   ///< BetaH charge
  int   beta_nhit;                ///< BetaH number of hits in charge computation
  float beta_q_err;               ///< BetaH charge error
  float beta_q_lay[4];            ///< BetaH layer charge
  float trd_chisq;                ///< TrdTrack Chi2/NDF of a linear fit
  int   trd_nseg;                 ///< TrdTrack number of segments in TrdTrack
  float trd_q;                    ///< TrdTrack charge
  bool  trdk_is_align_ok;         ///< TrdK alignment is loaded
  bool  trdk_is_calib_ok;         ///< TrdK calibration is loaded
  float trdk_q[3];                ///< TrdK charge (dE/dx&&delta-rays,dE/dx,delta-rays)
  float trdk_q_err[3];            ///< TrdK charge error
  int   trdk_like_valid;          ///< TrdK valid likelihood calculation
  float trdk_like_e;              ///< TrdK -loglikelihood for electron
  float trdk_like_p;              ///< TrdK -loglikelihood for proton
  float trdk_like_he;             ///< TrdK -loglikelihood for helium
  int   trdk_like_nhit;           ///< TrdK number of hits
  int   trdk_like_nhit_off;       ///< TrdK number of hits off the likelihood calculation
  float trdk_like_ampl_off;       ///< TrdK total amplitude of hits off the likelihood calculation
  float exthit_int[2][3];         ///< External hits extrapolation (TrdTrack or, if not found, BetaH) (L1,L9|x,y,z) [cm]
  float exthit_closest_coo[2][3]; ///< Coordinates of the external hits closest to the TrdTrack (or, if not found, closest to BetaH) (L1,L9|x,y,z) [cm]
  float exthit_closest_q[3][2];   ///< Charge of the external hits closest to the TrdTrack (or, if not found, closest to BetaH) (Old,Y.Jia,H.Liu|L1,L9) [cm]
  int   exthit_closest_status[2]; ///< Status of the external hits closest to the TrdTrack (or, if not found, closest to BetaH) (L1,L9)
  float exthit_largest_coo[2][3]; ///< Coordinates of the external hits with the largest charge (L1,L9|x,y,z) [cm]
  float exthit_largest_q[3][2];   ///< Charge of the external hits with the largest charge (Old,Y.Jia,H.Liu|L1,L9) [cm]
  int   exthit_largest_status[2]; ///< Status of the external hits with the largest charge (L1,L9)

  //! Calculate the average charge combining several TOF layers
  double GetSAQTof(int &n, double& rms, int pattern);

  NtpStandAlone(){}
  virtual ~NtpStandAlone(){}
  ClassDef(NtpStandAlone,3);
};

/** \class Event
All the informations related to a single event.
*/
class Event {

 public:

  NtpSHeader        *SHeader;
  NtpHeader         *Header;   ///< Header (run, event, time, prescaling, recon. counters, daq, trigger, orbit)
  NtpTrd            *Trd;      ///< TRD related info
  NtpTof            *Tof;      ///< TOF related info
  NtpTracker        *Tracker;  ///< Tracker related info
  NtpRich           *Rich;     ///< RICH related info
  NtpEcal           *Ecal;     ///< ECAL related info
  NtpAnti           *Anti;     ///< ACC related info
  NtpStandAlone     *SA;       ///< Stand-alone reconstruction related info
  NtpMCHeader       *MCHeader; ///< Monte Carlo truth related info

  ///< RTI related info
  RTIInfo           *RTI; //!

  ///< RICH occupancy calculation
  RichOccupancy     *RichOcc;

  ///< RichBDT calculator
  RichBDTMgr        *RichBDT; //!

  ///< TrackerBDT calculator
  TrackerBDTMgr     *TrackerBDT; //!

  ///< Classifier interface;
  ClassifierManager *Classifier; //!

  //! F. Dimiccoli RICH estimator
  Double_t GetRichBDT(bool debug=false);
  //! V. Formato Tracker estimator
  Double_t GetTrackerBDT(bool debug=false);
  //! Several other estimators (see Classifier.h for details)
  Double_t GetClassifier(int classifier_type);
  //! Get -LogLikelihood (see Classifier.h for details)
  Double_t GetMinusLogLikelihood(int ll_type);

  Event();
  virtual ~Event();
  ClassDef(Event,1);
};

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

  NtpCompact(){}
  virtual ~NtpCompact(){}
  bool IsStandalone();
  bool IsAnalysis();
  int IsInsideRich();

  int nParticle()    { return status%10; }
  int nAntiCluster() { return int(status/10)%10; } 
  int nBetaH()       { return int(status/100)%10; }
  int nTrTrack()     { return int(status/1000)%10; }
  int nTrRecHit()    { return int(status/10000)%100; }
  int nTrdCluster()  { return int(status/1000000)%100; } 
  int nTofClusterH() { return int(status/100000000)%100; } 

  ClassDef(NtpCompact,1);
};

#endif
