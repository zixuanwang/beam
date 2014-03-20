#include "Tracker.h"

namespace Beam{
	Tracker::Tracker(double tau_threshold_up, double tau_threshold_down, double level, double time) : m_tau_threshold_up(tau_threshold_up), m_tau_threshold_down(tau_threshold_down), m_level(level), m_time(time), m_signal(0){
	}

	Tracker::~Tracker(){

	}

	double Tracker::nextLevel(double time, double level){
		if (level < m_level){
			m_level += (time - m_time) / m_tau_threshold_down * (level - m_level);
			if (m_level < level){
				m_level = level;
			}
		}
		else{
			m_level += (time - m_time) / m_tau_threshold_up * (level - m_level);
			if (m_level > level){
				m_level = level;
			}
		}
		m_time = time;
		if (level > 3.5 * m_level){
			m_signal = 1;
		}
		if (level < 1.5 * m_level){
			m_signal = 0;
		}
		return m_level;
	}

	double Tracker::getLevel(){
		return m_level;
	}

	void Tracker::setLevel(double level){
		m_level = level;
	}

	int Tracker::classify(double time, double level){
		nextLevel(time, level);
		return m_signal;
	}
}
