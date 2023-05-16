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
