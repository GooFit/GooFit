#if !defined(FAKE_TH1F_H)
#define FAKE_TH1F_H 1

template<typename T>
class TH1 {
public:
   TH1(string const& name, string const& title, int numbins, int lowerlimit, int upperlimit) :
      name_(name), title_(title), numbins_(numbins), lowerlimit_(lowerlimit), upperlimit_(upperlimit),
      hist_(numbins + 2), stats_(false)
   {
      scale_ = (numbins_ - 1)/(upperlimit_ - lowerlimit_);
   }
   void SetStats(bool const stats) { stats_ = stats; }

   void Fill(T const val, T const weight=1.0) {
      int ival = (val - lowerlimit_)*scale_;
      if (ival < 0) {
	 ival = 0;			// off the left
      } else if (ival >= numbins_) {
	 ival = numbins_ + 1;		// off the right
      } else {
	 ++ival;			// we reserve [0] for underflow
      }
      hist_[ival] += weight;
   }

   void SetBinContent(int const i, T const val) {
      hist_[i] = val;
   }

   T GetBinContent(int const i) const {
      return hist_[i];
   }

   T GetBinCenter(int const i) const {
      return lowerlimit_ + (i - 1)/scale_;
   }

   int GetNumbins() const {
      return numbins_;
   }
private:
   std::string name_, title_;
   int numbins_, lowerlimit_, upperlimit_;
   std::vector<T> hist_;
   bool stats_;
   float scale_;
};

typedef TH1<float> TH1F;
#endif
