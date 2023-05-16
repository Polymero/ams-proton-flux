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
