type PositionEntry = { specCount: number; nonspecCount: number; strain: string };

type KmerHit = { entry: PositionEntry; spec: boolean };

export type ParsedReference = {
  kmers: Map<string, KmerHit>;
  positions: Map<string, PositionEntry>;
  strains: string[];
};

export function parseReference(tsv: string): ParsedReference {
  const kmers = new Map<string, KmerHit>();
  const positions = new Map<string, PositionEntry>();
  const strainSet = new Set<string>();

  const lines = tsv.split('\n');
  for (let i = 1; i < lines.length; i++) {
    const line = lines[i].trim();
    if (!line) continue;

    const [chromosome, position, specKmer, nonspecKmer, strain] = line.split('\t');
    const positionId = `${chromosome}:${position}`;

    let entry = positions.get(positionId);
    if (!entry) {
      entry = { specCount: 0, nonspecCount: 0, strain };
      positions.set(positionId, entry);
    }

    kmers.set(specKmer, { entry, spec: true });
    kmers.set(nonspecKmer, { entry, spec: false });
    strainSet.add(strain);
  }

  return { kmers, positions, strains: [...strainSet] };
}

export async function loadReference(url: string | URL): Promise<ParsedReference> {
  const response = await fetch(url);
  if (!response.ok) throw new Error(`Failed to fetch reference: ${response.status} ${url}`);
  return parseReference(await response.text());
}
