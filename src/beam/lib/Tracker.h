#ifndef TRACKER_H_
#define TRACKER_H_

namespace Beam{
	class Tracker {
	public:
		Tracker(double tau_threshold_up, double tau_threshold_down, double level, double time);
		~Tracker();
		double nextLevel(double time, double level);
		double getLevel();
		void setLevel(double level);
		int classify(double time, double level);
	private:
		double m_tau_threshold_up;
		double m_tau_threshold_down;
		double m_level;
		double m_time;
		int m_signal;
	};
}

#endif /* TRACKER_H_ */
