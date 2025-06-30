#include <cstdio>
#include <random>
#include <vector>

int main(int argc, char *argv[]) {
  if (argc != 3) {
    fprintf(stderr, "Usage: %s <number_of_gates> <output_file>\n", argv[0]);
    return 1;
  }

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> dis(0, 4);
  size_t N = atoll(argv[1]);
  std::vector<char> Gates(N);
  for (size_t i = 0; i < N; ++i) {
    int gate = dis(gen);
    Gates[i] = "HXYZS"[gate];
  }

  FILE *File = fopen(argv[2], "wb");
  if (!File) {
    perror("Failed to open file");
    return 1;
  }

  fwrite(&N, sizeof(size_t), 1, File);
  fwrite(Gates.data(), sizeof(char), N, File);
  fclose(File);

  return 0;
}
