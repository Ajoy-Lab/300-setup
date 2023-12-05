function new_phase = phase_wrap_pi_to_m_pi(phase)
    new_phase = phase;
    if phase > pi
        new_phase = phase - 2*pi;
    elseif phase < -pi
        new_phase = phase + 2*pi;
    else
    end
end