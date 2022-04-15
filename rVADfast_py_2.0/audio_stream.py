import os
import sys
import code
import pyaudio
import wave

# record stream audio, save it and then apply rVADfast to it
# usage: python3 audio_stream.py

CHUNK = 1024
FORMAT = pyaudio.paInt16
CHANNELS = 1
RATE = 44100

def record():

	p = pyaudio.PyAudio()

	stream = p.open(format=FORMAT,channels=CHANNELS,rate=RATE,input=True,frames_per_buffer=CHUNK)

	print("Start recording")
        print("2. Press Ctrl+C to stop the recording"
        print("3. rVAD will start")
        print("==================================================================\n")
        frames = []

	try:
		while True:
			data = stream.read(CHUNK)
			frames.append(data)

	except KeyboardInterrupt:
            print("Done recording: stored --> output.wav")

	except Exception as e:
		print(str(e))

	sample_width = p.get_sample_size(FORMAT)
	
	stream.stop_stream()
	stream.close()
	p.terminate()
	
	return sample_width, frames	

def record_to_file(file_path):
	wf = wave.open(file_path, 'wb')
	wf.setnchannels(CHANNELS)

	sample_width, frames = record()

	wf.setsampwidth(sample_width)
	wf.setframerate(RATE)
	wf.writeframes(b''.join(frames))
	wf.close()

if __name__ == '__main__':
        record_to_file('output.wav')
        print("rVAD running...")
        os.system("python3 rVAD_fast.py output.wav output.txt")
        print("Result written")
