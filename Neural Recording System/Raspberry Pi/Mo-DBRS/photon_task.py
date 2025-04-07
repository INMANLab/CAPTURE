
#!/usr/bin/python
import subprocess
import socket
import time
from rns import ResearchController
import sys
import signal
import RPi.GPIO as GPIO
import multiprocessing
import os
import datetime
import get_ntp 
import requests
import pyaudio
from rpi_ws281x import *

# LED strip configuration:
LED_COUNT      = 20      # Number of LED pixels.
LED_PIN        = 18      # GPIO pin connected to the pixels (must support PWM!).
LED_FREQ_HZ    = 800000  # LED signal frequency in hertz (usually 800khz)
LED_DMA        = 10      # DMA channel to use for generating signal (try 10)
LED_BRIGHTNESS = 255     # Set to 0 for darkest and 255 for brightest
LED_INVERT     = False   # True to invert the signal (when using NPN transistor level shift)
LED_CHANNEL    = 0
#LED_STRIP      = ws.SK6812_STRIP_RGBW
LED_STRIP      = ws.SK6812W_STRIP



# Define functions which animate LEDs in various ways.
def colorWipeLocal(strip, color, wait_ms=50):
    """Wipe color across display a pixel at a time."""
    for i in range(strip.numPixels()):
        strip.setPixelColor(i, color)
        strip.show()
        time.sleep(wait_ms/1000.0)


# NTP-Server for sync (first disable internal clock correction; temporary fix)
#os.system('sudo systemctl stop ntp.service')
goodclock = get_ntp.GoodClock()

# Command URL
#command_url = 'http://128.97.92.77'

# Define GPIOs
#
GPIO.setmode(GPIO.BOARD)
GPIO.setup(32, GPIO.IN, pull_up_down=GPIO.PUD_UP)
GPIO.setup(40, GPIO.OUT)
GPIO.setup(38, GPIO.OUT)
GPIO.output(40, GPIO.HIGH)
GPIO.output(38, GPIO.HIGH)

home_directory = sys.argv[1]



def StartStopBlink():

    n = 50 
    while n > 0:
        GPIO.output(40, True)
        GPIO.output(38, True)
        time.sleep(0.05)
        GPIO.output(40, False)
        GPIO.output(38, False)
        time.sleep(0.05)
        n = n - 1


# 50 ms pulse sent to GoPro LED.
#
def GoProLED(e, enable):
    while(1):
        e.wait()
        if enable.value:
            timestamp = goodclock.now()
            GPIO.output(40, GPIO.HIGH)
            time.sleep(0.05)
            GPIO.output(40, GPIO.LOW)
            TimestampQueue.put('GoPRO_VIDEO: '+str(timestamp) + '\n')
            print('GoPro LED Mark: ' + str(timestamp))
        else:
            e.clear()
            break
        e.clear()

# 50 ms pulse sent to Pupil LED.
#
def PupilLED(e, enable):
    while(1):
        e.wait()
        if enable.value:
            timestamp = goodclock.now()
            GPIO.output(38, GPIO.HIGH)
            time.sleep(0.05)
            GPIO.output(38, GPIO.LOW)
            TimestampQueue.put('PUPIL_VIDEO: '+str(timestamp) + '\n')
            print('Pupil LED Mark: ' + str(timestamp))
        else:
            e.clear()
            break
        e.clear()

# Sending RNS mark function.
#
def RNSMark(e,controller, msg, enable):
    while(1):
        e.wait()
        if enable.value:
            timestamp = goodclock.now()
            controller._send_message_check_reply_(msg, '')
            TimestampQueue.put('MARK: '+str(timestamp) + '\n')
            print('RNS Mark: ' + str(timestamp))
        else:
            e.clear()
            break
        e.clear()


n=0
stream=pyaudio.PyAudio().open(format=pyaudio.paInt8,channels=1,rate=44100,output=True)
square = chr(63)+ chr(0)

# Function sending an audio beep to GoPro audio input.
### Needs shorter audio file. Aplay delay? ###
#
def GoProAudio(e, enable):# controller, msg):
    while(1):
        e.wait()
        if enable.value:
            # Save timestamp
            timestamp = goodclock.now()
            for n in range(0,30,1): stream.write(square,50)
            TimestampQueue.put('GoPRO_AUDIO: '+str(timestamp) + '\n')
            print('Audio Mark: ' + str(timestamp))
        else:
            e.clear()
            break
        e.clear()

# Function writing timestamps from a queue to a file.
#
def TimestampWriter(dest_filename, some_queue, some_stop_token):
    with open(dest_filename, 'a+') as dest_file:
        while True:
            line = some_queue.get()
            if line == some_stop_token:
                return
            dest_file.write(line)



controller = 0
# If Research Accessories are connected open serial port for sending marks. Baud rate 9600 for old 300M device and 57600 for new 320 device.
#
try:
    controller = ResearchController(57600, 1)
except:
    print('\nCheck USB connection between RPi and Arduino!\n')
    print('\nContinuing without RNS!\n')
    controller = 0



# Define parallel processes
#
e = multiprocessing.Event()
taskEnabledFlag = multiprocessing.Value("i", 1)
msg = b't'

buttonEnabled = True
stopCounter = 0
STOP_TOKEN="STOP!!!"
TimestampQueue = multiprocessing.Queue()
rp_timestamp_str = str(datetime.datetime.now())
rp_timestamp_str = rp_timestamp_str.replace(' ', '_')
rp_timestamp_str = rp_timestamp_str.replace(':', '-')
rp_timestamp_str = rp_timestamp_str.replace('.', '_')

TimestampWriteProcess = multiprocessing.Process(target = TimestampWriter, args=(home_directory + "rp_timestamps/RP_marks_" + rp_timestamp_str + ".txt", TimestampQueue, STOP_TOKEN))
TimestampWriteProcess.start()



buttonEnabled = False
# Start NTP Time stamping process:
goodclock.run()
while goodclock.now() == None:
    print("We don't have a fix yet...")
    time.sleep(6)
print("Got a fix: "+str(goodclock.now()))

StartStopBlink()
timestamp = goodclock.now()
TimestampQueue.put('START: '+str(timestamp) + '\n')

# Run Matrix Creator script
#os.system('/home/pi/Downloads/MatrixCreatorAudioRecorder-master/build/microphone_array/mic_record_file &')
MatrixCreatorProcess = subprocess.Popen(['/home/pi/Downloads/MatrixCreatorAudioRecorder-master/build/microphone_array/mic_record_file'])


if controller is not 0:
    print ("Controller", controller)
    RNSProcess = multiprocessing.Process(target=RNSMark, args=(e,controller,msg, taskEnabledFlag))#controller, msg))
    RNSProcess.start()


GoProLEDProcess = multiprocessing.Process(target=GoProLED, args=(e, taskEnabledFlag,))
GoProLEDProcess.start()
PupilLEDProcess = multiprocessing.Process(target=PupilLED, args=(e,taskEnabledFlag,))
PupilLEDProcess.start()
GoProAudioProcess = multiprocessing.Process(target=GoProAudio, args=(e, taskEnabledFlag,))
GoProAudioProcess.start()



def ClosingRoutine():
        global taskEnabledFlag
        global e

        taskEnabledFlag.value = not taskEnabledFlag.value
        timestamp = goodclock.now()
        TimestampQueue.put('STOP: '+str(timestamp) + '\n')
        print("Interrupt! Closing app...")


        GPIO.output(40, GPIO.LOW)
        GPIO.output(38, GPIO.LOW)

        StartStopBlink()
#        NeoFinish(strip)

        time.sleep(20)

        print('Stopping processes...')
        e.set()
        time.sleep(2)
        e.clear()

        if controller is not 0:
                controller._send_message_check_reply_(msg, '')
                timestamp = goodclock.now()
                TimestampQueue.put('Close Mark: '+str(timestamp) + '\n')
                print('Close Mark: ' + str(timestamp))
                time.sleep(5)
                RNSProcess.terminate()
                controller.close()

        stream.close()
        MatrixCreatorProcess.terminate()
        pyaudio.PyAudio().terminate()
        GoProLEDProcess.terminate()
        PupilLEDProcess.terminate()
        GoProAudioProcess.terminate()
        print(STOP_TOKEN)
        TimestampQueue.put(STOP_TOKEN)
        time.sleep(2)
        TimestampWriteProcess.terminate()

        
        timestamp = goodclock.now()
    
        directory = 'matrix_mic_data/matrix_mic_' + str(timestamp) + '/'
        directory = directory.replace(' ', '_')
        directory = directory.replace(':', '-')
        directory = directory.replace('.', '_')
        directory = home_directory + directory
        print(directory)
        os.system('mkdir ' + directory)
        os.system('sudo mv ./mic_* ' + directory)


        GPIO.cleanup()
        sys.exit()


def buttonEvent_Falling(enabled):

    global stopCounter

    print('ButtonPress')

    if enabled:
        state = GPIO.input(40)
        if state:
            GPIO.output(40, False)
            GPIO.output(38, False)
            time.sleep(0.2)

        GPIO.output(40, True)
        GPIO.output(38, True)
        time.sleep(0.2)
        GPIO.output(40, False)
        GPIO.output(38, False)

    else:
        stopCounter = stopCounter + 1
        if(stopCounter > 5):
            stopCounter = 0
            ClosingRoutine()
            

GPIO.add_event_detect(32, GPIO.FALLING, bouncetime=300, callback=lambda x: buttonEvent_Falling(buttonEnabled))

# Ctrl-C for closing the script
#
def signal_handler(sig, frame):
    pass
#    ClosingRoutine()

signal.signal(signal.SIGINT, signal_handler)



# Create NeoPixel object with appropriate configuration.
#strip = Adafruit_NeoPixel(LED_COUNT, LED_PIN, LED_FREQ_HZ, LED_DMA, LED_INVERT, LED_BRIGHTNESS, LED_CHANNEL, LED_STRIP)
# Intialize the library (must be called once before other functions).
#strip.begin()


'''
response = "false"
while response == "false":
    try:
        r = requests.get(command_url)
        response = r.content
    except requests.exceptions.RequestException as exc:
        print exc
    print("Response was false..")
    time.sleep(5)
'''
sync_period = 60
rns_period = 240
pre_ssr_period = 2
ssr_counter = 0
if controller is not 0:
    #controller.send_ssr()
    time.sleep(5)

try:
    while(1):
        # Trigger sync pulses
        if taskEnabledFlag.value:
            e.set()

            # Save timestamp
            timestamp = goodclock.now()
            TimestampQueue.put('GLOBAL: '+str(timestamp) + '\n')
            print('GLOBAL: ' + str(timestamp))
            #print(datetime.datetime.now())
            if taskEnabledFlag.value:
                #colorWipeLocal(strip, Color(255, 255, 255, 255), 0)
                #time.sleep(0.05)
                #colorWipeLocal(strip, Color(0, 0, 0, 0), 0)
                #ssr_counter = ssr_counter + 1
                #if (controller is not 0):
                #    if (ssr_counter + 1) * sync_period >= rns_period:
                        #time.sleep(pre_ssr_period)
                        #controller.send_ssr()
                        #timestamp = goodclock.now()
                        #TimestampQueue.put('SSR: '+str(timestamp) + '\n')
                        #print('SSR: ' + str(timestamp))
                        #ssr_counter = 0

                time.sleep(sync_period)
        else:
            break
        
    GoProLEDProcess.join()
    PupilLEDProcess.join()
    GoProAudioProcess.join()
except:
    print('EXCEPTION')
    #ClosingRoutine()
