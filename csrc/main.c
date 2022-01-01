#include <stdio.h>
#include <stdint.h>
#include "main.h"

#define BUFF_SIZE 4096

int main() {
    votrax_start();

    int16_t soundBuffer[BUFF_SIZE*2];

    votraxsc01_w(28); // G

    Votrax_Update(soundBuffer, BUFF_SIZE-1);

    votraxsc01_w(33); // AY

    return 1;
}

