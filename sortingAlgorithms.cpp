#include <iostream>
#include <fstream> //used for writing to files
#include <bits/stdc++.h> //used for readnig files
#include <chrono> //used for timing

using namespace std;

void printArray(int arr[], int size){
    for (int i = 0; i<size; i++)
        cout << arr[i] << ' ';
    cout << endl;
}

void selectionSort(int a[], int size){
  for (int i=0; i<size-1; i++){
    int minKey = i;
    // get minKey of suffix
    for (int j=i+1; j<size; j++){
      if (a[j] < a[minKey])
        minKey=j;
    }

    // swap, if needed
    if (i != minKey){
      a[i]  = a[i] + a[minKey];
      a[minKey] = a[i] - a[minKey];
      a[i] = a[i] - a[minKey];
    }
  }
}

void insertionSort(int a[], int size){
  for (int i = 1; i < size; ++i) {
    int key = a[i];
    int j = i - 1;

    while (j >= 0 && a[j] > key) {
      a[j + 1] = a[j];
      j--;
    }
    a[j + 1] = key;
  }
}

void merge(int a[], int left, int mid, int right){
  int n1 = mid - left + 1;
  int n2 = right - mid;
  int L[n1];
  int R[n2];
  //make auxilliary arrays
  for (int i = 0; i < n1; i++)
    L[i] = a[left + i];
  for (int j = 0; j < n2; j++)
    R[j] = a[mid + 1 + j];

  int i = 0;
  int j = 0;
  int k = left;

  while (i < n1 && j < n2) {
    if (L[i] <= R[j]) {
      a[k] = L[i];
      i++;
    }
    else {
      a[k] = R[j];
      j++;
    }
    k++;
  }
  while (i < n1) {
    a[k] = L[i];
    i++;
    k++;
  }
  while (j < n2) {
    a[k] = R[j];
    j++;
    k++;
  }
}

void mergeSort(int a[], int left, int right){
  if (left >= right)
    return;
  int mid = left + (right - left) / 2; //since datatype is int, drops floating point digits
  mergeSort(a, left, mid);
  mergeSort(a, mid + 1, right);
  merge(a, left, mid, right);
}

int partition(int a[], int left, int right){
  int pivot = a[right];
  int i = left-1;
  for(int j = left; j < right; j++){
    if (a[j] <= pivot){
      i++;
      int temp = a[i];
      a[i] = a[j];
      a[j] = temp;
    }
  }
  i++;
  int temp = a[i];
  a[i]=a[right];
  a[right]=temp;
  return i;
}

void quickSort(int a[], int left, int right){
  if(left >= right)
    return;
  int p = partition(a, left, right);
  quickSort(a, left, p-1);
  quickSort(a, p+1, right);
}

int main() {
  int num_iter = 20;
  int max_size = 1000;
  ofstream myFile;
  myFile.open("quickSortTimes.txt");
  for(int size = 1; size <= max_size; size++){ //loop over input sizes
    double t=0;
    for (int iter = 1; iter <= num_iter; iter++){ //different inputs for each size
      //Declaring Variables
      int count = 0;
      int num;
      int arr[size];

      //Reading Input Array
      ifstream fp("Permutations/permutation"+to_string(iter)+".txt");
      if (!fp.is_open()) { // check if the file is opened
        cerr << "error opening the file";
        return 1;
      }
      while(fp >> num) {
        arr[count++] = num;
        if (count == size)
          break;
      }
      fp.close();
      
      //Sort and add time
      auto start = std::chrono::high_resolution_clock::now();
      quickSort(arr, 0, size-1);
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> duration = end - start;
      t = t + duration.count();
    }
    //print average time to file
    myFile << size << "," << t/num_iter << endl;
  }
  myFile.close();
  cout << "Program Ended";
  return 0;
}