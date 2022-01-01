#include <stdio.h>
#include <stdint.h>
#include "main.h"

#include <pulse/simple.h>

#define BUFF_SIZE 4096

int main() {
    pa_sample_spec ss;

    ss.format = PA_SAMPLE_S16NE;
    ss.channels = 1;
    ss.rate = 9000;

    pa_simple* s = pa_simple_new(NULL,               // Use the default server.
                      "Votrax Testing",           // Our application's name.
                      PA_STREAM_PLAYBACK,
                      NULL,               // Use the default device.
                      "Music",            // Description of our stream.
                      &ss,                // Our sample format.
                      NULL,               // Use default channel map
                      NULL,               // Use default buffering attributes.
                      NULL               // Ignore error code.
    );



    votrax_start();

    int16_t soundBuffer[BUFF_SIZE*2];

    votraxsc01_w(28); // G

    Votrax_Update(soundBuffer, BUFF_SIZE-1);

    //for(int x=0; x<BUFF_SIZE; x++) {
    //    printf("%d: %x\n",x, soundBuffer[x]);
    //}

    int error;
    int res;

    res = pa_simple_write(s, soundBuffer, BUFF_SIZE, &error);
    printf("WRITE ERROR %d VAL %d\n", res, error);

    res = pa_simple_drain(s, &error);
    printf("DRAIN ERROR %d VAL %d\n", res, error);

    votraxsc01_w(33); // AY
    Votrax_Update(soundBuffer, BUFF_SIZE-1);

    res = pa_simple_write(s, soundBuffer, BUFF_SIZE, &error);
    printf("WRITE ERROR %d VAL %d\n", res, error);

    res = pa_simple_drain(s, &error);
    printf("DRAIN ERROR %d VAL %d\n", res, error);
    
    pa_simple_free(s);

    return 1;
}

