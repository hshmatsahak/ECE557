#include "Timer.h"
#include <BasicLinearAlgebra.h>
using namespace BLA;

/** -----------------------------------------------------------
 * SECTION 1
 * 
 * @students
 * TODO: Add your system matrices below
 *
* ----------------------------------------------------------- **/
Matrix<2, 1> x_hat; //observer states, initialized to 0 in setup
Matrix<2, 1> x_hat_dot;
Matrix<1, 1> u = {0}; //Initial control input value
double y; // measured output

/** YOUR MATRICES HERE 
 *  Below is an example of how to declare a 2 by 2 identity matrix,
 *  you can delete this identityMatrix as it is not used in the
 *  rest of the code
*/
//Matrix<2, 2> MatrixA = {
//-107.212276214834,  1,
//-3828.99902538576, -80
//};
//
//Matrix<2, 3> MatrixB = {
//107.212276214834,  0, 0,
//2228.99902538576,  1600,  67.2122762148338
//};
//
//Matrix<1, 2> MatrixC = {
//-760.938314509123,  -31.9652476107732
//};
//
//Matrix<1, 3> MatrixD = {
//0,  760.938314509123,  31.9652476107732
//};

Matrix<2, 2> MatrixA = {
-87.2122762148338,  1,
-3884.75350108908, -100
};

Matrix<2, 3> MatrixB = {
87.2122762148338,  0, 0,
1384.75350108908,  2500,  87.2122762148338
};

Matrix<1, 2> MatrixC = {
-1188.96611642050,  -41.4769765421373
};

Matrix<1, 3> MatrixD = {
0,  1188.96611642050,  41.4769765421373
};

Matrix<2, 2> identityMatrix = {
  1, 0,
  0, 1
};

//Signal Generator parameters
double f_sg = 0.1; // signal generator frequency (Hz)
double T_sg = 1/f_sg; //signal generator period

double a_sg = 0.03; //signal generator amplitude (m)
double z_d = -a_sg; // reference signal initial value

/** -----------------------------------------------------------
 * SECTION 2
 * Do NOT modify
 * 
 * Sets up timer, encoder, and motor pin constants
 * ----------------------------------------------------------- **/
Timer t; 
double T_sample = 0.001;  // Controller sample time in seconds
double T_plot = 0.03;     // writing to serial every T_plot seconds
char serial_delim = ',';  // delimeter for serial interface communications

// cart position encoder
enum PinAssignments{ // cart encoder pins
  encoderPinA = 2,  // brown wire, blk wire gnd 
  encoderPinB = 3,  // green wire, red wire 5v
};
volatile int encoderPos = 0; //count cart ticks
boolean A_set = false; //Boolean variables for determining encoder position changes
boolean B_set = false;
const float pi = 3.14;
double K_encoder = 2.2749 * 0.00001; // constant to map encoder clicks to position

// motor control
int PWM_PIN = 11;   // pin to send the magnitude of applied motor voltage
int DIR_PIN = 8;    // pin to send the polarity of applied motor voltage
int duty_cycle;   // duty cycle of pwm to send to motor

/** -----------------------------------------------------------
 * SECTION 3
 * Do NOT modify
 * 
 * Setup function which sets the mode for Arduino pins
 * ----------------------------------------------------------- **/
void setup() {
  Serial.begin (115200); // open serial connection for writing
  x_hat.Fill(0); // initialize observer states to 0
  
  configurePins(); // sets up Arduino pin modes
  
  t.every(T_sample*1000, outputFeedbackControl);  // update observer and controller every T_sample seconds 
  t.every(T_plot*1000, writeToSerial);            // perform writeToSerial every T_plot seconds
  t.every(T_sg*1000/2, toggleSquareWave);         // switches signal amplitude between + and -
                                                  // t.every takes in ms
}

/** -----------------------------------------------------------
 * SECTION 4
 * Do NOT modify
 * 
 * main loop() which executes over and over again on the Arduino
 * ----------------------------------------------------------- **/
void loop(){
  t.update(); // update timer, which triggers scheduled functions 
              // according the t.every(...) set in the setup() function
}

/** -----------------------------------------------------------
 * SECTION 5
 * 
 * @students
 * TODO: Add the steps to update the observer state estimate 
 * and the feedback control input in the function 
 * outputFeedbackControl()
 *
 * ----------------------------------------------------------- **/
void outputFeedbackControl(){
    /* YOUR CODE HERE
     *  make sure to update the value of the global variable u 
     *  before the function call to mapControlInputToMotorShield()
    */

    y = K_encoder * encoderPos;
    
//    Matrix <2, 1> x_d = {
//      z_d,
//      0
//    };
    
    Matrix <3, 1> u_hat = {
      y,
      z_d,
      0
    };
    
    x_hat_dot = MatrixA * x_hat + MatrixB * u_hat;
    x_hat = x_hat + x_hat_dot * T_sample; // Euler discretization step
    
    u = MatrixC * x_hat + MatrixD * u_hat;
    
    mapControlInputToMotorShield();
}  

/** -----------------------------------------------------------
 * SECTION 6
 * Do NOT modify
 * 
 * Utility functions for creating the square wave reference signal, 
 * writing to the serial port, and interfacing with the motor shield
 * ----------------------------------------------------------- **/
void toggleSquareWave(){
  z_d = -z_d;
}
  
void writeToSerial(){ 
  Serial.print(y, 5); 
  Serial.print(serial_delim);
  Serial.print(z_d, 5);
  Serial.print(serial_delim);
  Serial.print(x_hat(0), 5);

  Serial.println();
}

void mapControlInputToMotorShield(){
  //Mapping between required voltage and 11.75V motor shield
  duty_cycle = round(u(0)/11.75*1024); 
  if (duty_cycle > 512 ){ //Saturation to not exceed motor voltage limits of 5.875 V
    duty_cycle=512; //motor moves left
  } else if (duty_cycle < -512 ){
    duty_cycle=-512;//motor moves right
  } 
  
  if (duty_cycle > 0){ // Update the direction of motor based on the sign of duty_cycle  
    digitalWrite(DIR_PIN, HIGH);
    OCR1A = duty_cycle;
  } else if (duty_cycle <= 0){
    digitalWrite(DIR_PIN, LOW);
    OCR1A = -duty_cycle;
  }
}

void configurePins(){
  pinMode(encoderPinA, INPUT); // cart encoder pins
  pinMode(encoderPinB, INPUT); 
  
  digitalWrite(encoderPinA, HIGH);  // turn on pullup resistor
  digitalWrite(encoderPinB, HIGH);  // turn on pullup resistor
  attachInterrupt(0, doEncoderA, CHANGE);
  attachInterrupt(1, doEncoderB, CHANGE);

  pinMode(PWM_PIN, OUTPUT); // motor pins
  pinMode(DIR_PIN, OUTPUT); 
  
  //Clock setup to increase PWM frequency controlling applied voltage on PIN 11
  TCCR1A = _BV(COM1A1) | _BV(WGM21) | _BV(WGM20);
  TCCR1B = _BV(CS10);
  OCR1A = 0;
}

/** -----------------------------------------------------------
 * SECTION 7
 * Do NOT modify
 * 
 * Interrupts for reading cart position encoders and updating
 * encoderPos
 * ----------------------------------------------------------- **/
void doEncoderA(){ // Interrupt on A changing state
  A_set = digitalRead(encoderPinA) == HIGH; // Test transition
  encoderPos += (A_set != B_set) ? +1 : -1; // and adjust counter + if A leads B
}

void doEncoderB(){ // Interrupt on B changing state
  B_set = digitalRead(encoderPinB) == HIGH; // Test transition
  encoderPos += (A_set == B_set) ? +1 : -1; // and adjust counter + if B follows A 
}
