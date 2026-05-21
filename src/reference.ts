type PositionEntry = { specCount: number; nonspecCount: number; strain: string };

export type ParsedReference = {
  specKmers: Map<string, string>;
  nonspecKmers: Map<string, string>;
  positions: Map<string, PositionEntry>;
  strains: string[];
};

export function parseReference(tsv: string): ParsedReference {
  const specKmers = new Map<string, string>();
  const nonspecKmers = new Map<string, string>();
  const positions = new Map<string, PositionEntry>();
  const strainSet = new Set<string>();

  const lines = tsv.split('\n');
  for (let i = 1; i < lines.length; i++) {
    const line = lines[i].trim();
    if (!line) continue;

    const [chromosome, position, specKmer, nonspecKmer, strain] = line.split('\t');
    const positionId = `${chromosome}:${position}`;

    specKmers.set(specKmer, positionId);
    nonspecKmers.set(nonspecKmer, positionId);
    strainSet.add(strain);

    if (!positions.has(positionId)) {
      positions.set(positionId, { specCount: 0, nonspecCount: 0, strain });
    }
  }

  return { specKmers, nonspecKmers, positions, strains: [...strainSet] };
}

export async function loadReference(url: string | URL): Promise<ParsedReference> {
  const response = await fetch(url);
  if (!response.ok) throw new Error(`Failed to fetch reference: ${response.status} ${url}`);
  return parseReference(await response.text());
}
