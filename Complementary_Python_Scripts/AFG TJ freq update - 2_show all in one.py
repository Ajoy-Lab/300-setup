import pyvisa
import time
import socket
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

ch= 1
volts = 0.15

silent = True


testmode = False
def main():
    if not testmode:
        #Replace this address:
        rm = pyvisa.ResourceManager()
        list_available_devices(rm)
        
        visa_address_RF = 'USB0::0x0699::0x035A::B011535::INSTR'       # Replace with your instrument's address
        visa_address_TJ = "USB0::0x0699::0x0355::C019986::INSTR"
        visa_address_AC = 'USB0::0x0699::0x0355::C019987::INSTR'
        afgRF = TektronixAFG31000(visa_address_RF, rm)                       # Create an instance of the TektronixAFG31000 class
        afgTJ = TektronixAFG31000(visa_address_TJ, rm)                       # Create an instance of the TektronixAFG31000 class
        afgAC = TektronixAFG31000(visa_address_AC, rm)
        channel = 1

        afgAC.identify()
        afgRF.identify()                                              # Print the instrument identification
        afgTJ.identify()                                              # Print the instrument identification
        
        # afgAC.afg.write('*RST')
        # afgRF.reset()
        # input()
        
        # f_Larmor = 75380000
        # cmds = [
        #     'SOUR1:FUNC SIN',
        #     'SOUR2:FUNC SIN',
        #     f"SOUR1:FREQ {f_Larmor}",
        #     f"SOUR2:FREQ {f_Larmor+2*500}",
        #     'SOUR1:VOLT 0.25',
        #     'SOUR2:VOLT 0.25',
        #     'SOUR1:VOLT:OFFS 0',
        #     'SOUR2:VOLT:OFFS 0',
        #     'SOUR1:FM:STAT ON',
        #     'SOUR2:FM:STAT ON',
        #     'SOUR1:FM:INT:FUNC TRI',
        #     'SOUR2:FM:INT:FUNC TRI',
        #     'SOUR1:FM:INT:FREQ 2.5Hz',
        #     'SOUR2:FM:INT:FREQ 2.5Hz',
        #     'SOUR1:FM:DEV 500Hz',
        #     'SOUR2:FM:DEV 500Hz'
        # ]
        
        # for cmd in cmds:
        #     afgRF.afg.write(cmd)
            
        # input()
        
        # afgAC.afg.write('BURSt:STATE ON')
        # afgAC.afg.write('SOURce1:BURSt:IDLE DC')
        # afgAC.afg.write('SOURce1:BURSt:INFInite:REARm')
        # afgAC.afg.write('SOURce1:BURSt:MODE GATE')
        # resp = afgAC.afg.write('SOUR1:FUNC:EFIL M:/.tfwx')
        # afgAC.afg.write('SOUR1:FUNC:SHAP EFIL')
        # afgAC.afg.write('SOUR1:FREQ 0.223869232')
        # afgAC.afg.write('SOUR1:VOLT 0.25')
        # afgAC.afg.write('SOUR1:VOLT:OFFS 0')
        # afgAC.afg.write('SOUR1:BURSt:NCYCles 1')
        # afgAC.afg.write('SOURce1:BURSt:INFInite:REARm')
        # input(f"response: {resp}")
        #afgRF.set_volts(channel=ch, voltage=volts)
        #afg.configure(frequency='1.552838kHz', voltage=0.5, channel=channel, phase=0)     # Configure the function generator

    else:
        afgRF = None
        afgTJ = None
    # while True:
    # try: 
    while True:
            response, code = wait_for_udp_packet(host = 'localhost', port = 12345, expected_message = 'start_output')
            if code == 2:
                rf_freq = float(response)
                afgRF.set_freq(rf_freq)
                print(f"RF freq set to {rf_freq}")
                continue
            frequency = float(response)
            time.sleep(1.25) #-0.07
            afgTJ.afg.write(f'SOUR{channel}:FREQ {frequency-2}') # use this for orbits
            time.sleep(0.75-0.002)
            afgTJ.afg.write(f'SOUR{channel}:FREQ {frequency}') 
            time.sleep(2)
            # afgAC.start_output()
            afgAC.afg.write(f'SOUR{channel}:FUNC SQU')
            time.sleep(1)
            afgRF.start_output()
            time.sleep(1)
            afgAC.set_volts(voltage = 0.001)
            afgAC.set_freq(0.01)
            
            # # afgRF.start_output(channel=2)
            # afgRF.start_output()
            # time.sleep(0.5)
            # # afgRF.start_output(channel=2)
            # # time.sleep(3.5)
            time.sleep(30)
            print(f"All set. Response was: {response} Hz")
            afgTJ.stop_output()
            afgAC.stop_output()
            afgRF.stop_output()
            afgRF.stop_output(channel=2)
    # except:
    #     afgTJ.stop_output()               # Turn off the output when the script is interrupted
    #     afgAC.stop_output()
    #     afgRF.stop_output()
    #     print("Outputs turned off.")
    #     afgTJ.close()  
    #     afgRF.close() 
    
    time.sleep(1)
    if not testmode:
        afgTJ.close()                         # Close the connection
        afgAC.close()
        afgRF.close() 



class TektronixAFG31000:
    def __init__(self, visa_address, rm):
        self.rm = rm
        self.afg = self.rm.open_resource(visa_address)

    def identify(self):
        """Print the instrument identification."""
        print(self.afg.query('*IDN?'))

    def configure(self, channel=1, waveform='DC', frequency='1kHz', voltage=1.0, phase = 0):
        """Configure the function generator with initial settings."""
        self.afg.write(f'SOUR{channel}:FUNC {waveform}')  # Set function to specified waveform
        #self.afg.write(f'SOUR{channel}:FREQ {frequency}') # Set frequency to specified value
        self.afg.write(f'SOUR{channel}:VOLT {voltage}')   # Set initial voltage amplitude
        #self.afg.write(f'SOUR{channel}:PHAS {phase}DEG')

    def set_freq(self, frequency, channel = 1):
        self.afg.write(f'SOUR{channel}:FREQ {frequency}') # Set frequency to specified value

    def set_offset(self, channel = 1, voltage = 1.0):
        #self.afg.write(f'SOUR{channel}:VOLT {voltage}')   # Set voltage amplitude
        self.afg.write(f'SOUR1:VOLT:LEV:IMM:OFFS {voltage}V')
    
    def set_volts(self, channel = 1, voltage = 1.0):
        self.afg.write(f'SOUR{channel}:VOLT {voltage}')   # Set voltage amplitude

    def get_volts(self, channel = 1):
        return float(self.afg.query(f'SOUR{channel}:VOLT?'))

    def set_phase_offset(self, channel = 1, phase = 0):
        None#self.afg.write(f'SOUR{channel}:PHAS {phase}DEG')

    def start_output(self, channel=1):
        """Start the output of the function generator."""
        self.afg.write(f'OUTP{channel} ON')

    def stop_output(self, channel=1):
        """Stop the output of the function generator."""
        self.afg.write(f'OUTP{channel} OFF')

    def close(self):
        """Close the connection to the function generator."""
        self.afg.close()
    
    def reset(self):
        self.afg.write("*RST")

def wait_for_udp_packet(host, port, expected_message):
    code = None
    """Wait for a specific UDP packet before starting the output."""
    sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    sock.bind((host, port))
    print(f"Waiting for UDP packet on {host}:{port}...")

    while True:
        data, addr = sock.recvfrom(1024)  # Buffer size is 1024 bytes
        if not silent: print(f"Received message: {data} from {addr}")
        data_dec = data.decode()
        if data_dec[:len(expected_message)] == expected_message:
            if not silent: print("Expected UDP packet received. Starting output.")
            code = 1
            break
        elif data_dec[:len(expected_message)] == 'set_rf_freq_':
            if not silent: print("Expected UDP packet received. Starting output.")
            code = 2
            break

    sock.close()
    assert(code is not None)
    try:
        return data_dec[len(expected_message):], code
    except:
        return None, None

def list_available_devices(rm):
    """List all available devices in the resource manager."""
    resources = rm.list_resources()
    print("Available devices:")
    for resource in resources:
        print(resource)

if __name__ == "__main__":
    main()





