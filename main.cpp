#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <string>
#include <stdexcept>
#include <cctype>
#include <limits>

enum class StrandLogic {
	CodingStrand,
	TemplateStrand
};

const std::unordered_map<char, char> TemplateStrandMap = {
	{'A', 'U'}, {'T', 'A'}, {'C', 'G'}, {'G', 'C'},
	{'R', 'Y'}, {'Y', 'R'}, {'S', 'S'}, {'W', 'W'},
	{'K', 'M'}, {'M', 'K'}, {'B', 'V'}, {'D', 'H'},
	{'H', 'D'}, {'V', 'B'}, {'N', 'N'}
};

const std::unordered_map<char, std::vector<char>> AmbiguousBasesRNA = {
	{'A', {'A'}}, {'U', {'U'}}, {'C', {'C'}}, {'G', {'G'}},
	{'R', {'A', 'G'}}, {'Y', {'C', 'U'}}, {'S', {'G', 'C'}}, {'W', {'A', 'U'}},
	{'K', {'G', 'U'}}, {'M', {'A', 'C'}}, {'B', {'C', 'G', 'U'}},
	{'D', {'A', 'G', 'U'}}, {'H', {'A', 'C', 'U'}}, {'V', {'A', 'C', 'G'}},
	{'N', {'A', 'U', 'C', 'G'}}
};

// generates all possible RNA sequences (recursively)
void generateAllRNA(const std::string& rna, std::string& current, size_t index, std::vector<std::string>& results) {
	if (index == rna.size()) {
		results.push_back(current);
		return;
	}
	char base = rna[index];
	auto it = AmbiguousBasesRNA.find(base);
	if (it != AmbiguousBasesRNA.end()) {
		for (char rnaBase : it->second) {
			current.push_back(rnaBase);
			generateAllRNA(rna, current, index + 1, results);
			current.pop_back();
		}
	}
	else {
		current.push_back(base);
		generateAllRNA(rna, current, index + 1, results);
		current.pop_back();
	}
}

// wrapper function to start RNA generation
std::vector<std::string> getAllRNA(const std::string& rna) {
	std::vector<std::string> results;
	std::string current;
	current.reserve(rna.size());
	generateAllRNA(rna, current, 0, results);
	return results;
}

// Converts DNA to RNA based on strand logic (coding or template)
std::string dnaToRna(const std::string& dna, StrandLogic logic) {
	std::string rna;
	rna.reserve(dna.size());
	for (char base : dna) {
		if (logic == StrandLogic::CodingStrand) {
			rna += (base == 'T') ? 'U' : base;
		}
		else if (logic == StrandLogic::TemplateStrand) {
			auto it = TemplateStrandMap.find(base);
			if (it != TemplateStrandMap.end()) {
				rna += it->second;
			}
			else {
				rna += '?';
			}
		}
	}
	return rna;
}

// loads an entire text file into a string
std::string loadFile(const std::string& filePath) {
	std::ifstream file(filePath);
	if (!file) throw std::runtime_error("Failed to open file: " + filePath);
	std::ostringstream buffer;
	buffer << file.rdbuf();
	return buffer.str();
}

// extracts the desc line from the FASTA file
std::string getFastaDescription(const std::string& content) {
	std::istringstream stream(content);
	char firstChar;
	stream >> firstChar;
	if (firstChar != '>') {
		throw std::runtime_error("FASTA header should start with '>'");
	}
	std::string description;
	std::getline(stream, description);
	return description;
}

// Extracts the DNA sequence from the FASTA content (all lines after the header)
std::string extractDNA(const std::string& content) {
	std::istringstream stream(content);
	std::string line;
	std::getline(stream, line);
	std::string sequence;
	while (std::getline(stream, line)) {
		for (char c : line) {
			if (!std::isspace(static_cast<unsigned char>(c))) {
				sequence += std::toupper(c);
			}
		}
	}
	return sequence;
}

// ask user with a yes or no question (helper)
bool askYesNo(const std::string& question) {
	std::string answer;
	std::cout << question << " (y/n): ";
	std::getline(std::cin, answer);
	return !answer.empty() && (answer[0] == 'y' || answer[0] == 'Y');
}

// saves all RNA sequences to a file
void saveToFile(const std::vector<std::string>& rnaList, const std::string& fileName) {
	std::ofstream out(fileName);
	if (!out) {
		std::cerr << "Failed to write to file: " << fileName << "\n";
		return;
	}
	for (size_t i = 0; i < rnaList.size(); ++i) {
		out << "RNA [" << i + 1 << "]: " << rnaList[i] << "\n";
	}
	std::cout << "Results saved to: " << fileName << "\n";
}

int main() {
	try {
		std::string filePath;
		std::cout << "Enter path to FASTA file: ";
		std::getline(std::cin, filePath);

		std::string content = loadFile(filePath);
		std::string description = getFastaDescription(content);
		std::string dna = extractDNA(content);

		std::cout << "\n===== FASTA Description =====\n" << description << "\n";
		std::cout << "\n===== DNA Sequence =====\n" << dna << "\n";

		std::cout << "\nUse which strand for RNA transcription?\n";
		std::cout << "1. Coding Strand\n";
		std::cout << "2. Template Strand\n";
		std::cout << "Choice: ";
		std::string choice;
		std::getline(std::cin, choice);

		StrandLogic logic = (choice == "2") ? StrandLogic::TemplateStrand : StrandLogic::CodingStrand;
		std::string rnaBase = dnaToRna(dna, logic);

		std::vector<std::string> rnaVariants = getAllRNA(rnaBase);

		std::cout << "\n===== RNA Variants (First 10 or less) =====\n";
		size_t displayLimit = std::min<size_t>(10, rnaVariants.size());
		for (size_t i = 0; i < displayLimit; ++i) {
			std::cout << "RNA [" << i + 1 << "]: " << rnaVariants[i] << "\n";
		}

		if (rnaVariants.size() > 10) {
			std::cout << "...and " << (rnaVariants.size() - 10) << " more variants.\n";
			if (askYesNo("Show all RNA variants?")) {
				for (size_t i = 10; i < rnaVariants.size(); ++i) {
					std::cout << "RNA [" << i + 1 << "]: " << rnaVariants[i] << "\n";
				}
			}
		}

		if (askYesNo("Do you want to save the RNA sequences to a file?")) {
			std::string outputFile;
			std::cout << "Enter output file name: ";
			std::getline(std::cin, outputFile);
			saveToFile(rnaVariants, outputFile);
		}

		std::cout << "\n===== Summary =====\n";
		std::cout << "Total RNA sequences generated: " << rnaVariants.size() << "\n";
		std::cout << "Original DNA length: " << dna.length() << "\n";

	}
	catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}

	return 0;
}
