export interface EnergyUpgrade {
  id: string;
  name: string;
  subtitle: string;
  description: string;
  savings: string;
  costRange: string;
  payback: string;
  enovaSupport: string;
  color: string;
  highlightColor: string;
}

export const UPGRADES: Record<string, EnergyUpgrade> = {
  roof: {
    id: "roof",
    name: "Loftsisolasjon",
    subtitle: "Tak og loft",
    description:
      "God isolasjon i taket er blant de mest lønnsomme energitiltakene. Varmen stiger – stopp det med riktig isolasjon og spar på oppvarmingskostnadene.",
    savings: "Spar opptil 25 % energi",
    costRange: "20 000 – 80 000 kr",
    payback: "5–10 år",
    enovaSupport: "Opptil 10 000 kr",
    color: "#6B3A2A",
    highlightColor: "#D97706",
  },
  walls: {
    id: "walls",
    name: "Veggsisolasjon",
    subtitle: "Vegger og kledning",
    description:
      "Etterisolering av vegger reduserer varmetapet betraktelig. Kombinert med ny kledning gir det både energibesparelser og fornyet fasade.",
    savings: "Spar opptil 20 % energi",
    costRange: "50 000 – 250 000 kr",
    payback: "10–20 år",
    enovaSupport: "Opptil 25 000 kr",
    color: "#D6CFC4",
    highlightColor: "#F59E0B",
  },
  windows: {
    id: "windows",
    name: "Vindusskifte",
    subtitle: "Vinduer og glass",
    description:
      "Bytte til energivinduer med tre-lags glass reduserer varmetapet kraftig og øker komforten – uten kulderas og trekk.",
    savings: "Spar opptil 15 % energi",
    costRange: "3 000 – 8 000 kr per vindu",
    payback: "8–15 år",
    enovaSupport: "Opptil 1 000 kr per vindu",
    color: "#BAE6FD",
    highlightColor: "#3B82F6",
  },
  door: {
    id: "door",
    name: "Ytterdør",
    subtitle: "Dør og inngangsparti",
    description:
      "En isolert ytterdør hindrer trekkfulle inngangspartier og kuldebroer. Kombinert med ny karm og tettelister gir det merkbar effekt.",
    savings: "Komfort + energibesparelse",
    costRange: "8 000 – 25 000 kr",
    payback: "10–20 år",
    enovaSupport: "Inngår i heltiltaksstøtte",
    color: "#1E293B",
    highlightColor: "#818CF8",
  },
  foundation: {
    id: "foundation",
    name: "Gulvisolasjon",
    subtitle: "Gulv og grunnmur",
    description:
      "Isolasjon under gulv og i krypkjeller fjerner kalde gulv og reduserer energibehovet – gjerne kombinert med andre arbeider.",
    savings: "Spar opptil 10 % energi",
    costRange: "15 000 – 60 000 kr",
    payback: "7–14 år",
    enovaSupport: "Opptil 8 000 kr",
    color: "#94A3B8",
    highlightColor: "#10B981",
  },
  solar: {
    id: "solar",
    name: "Solcelleanlegg",
    subtitle: "Fornybar strømproduksjon",
    description:
      "Solceller på taket produserer strøm og reduserer strømregningen. Overskuddsstrøm selges tilbake til nettet. Enklere å installere enn du tror.",
    savings: "Produser 30–70 % av strømbehovet",
    costRange: "80 000 – 200 000 kr",
    payback: "8–12 år",
    enovaSupport: "Opptil 25 000 kr",
    color: "#1E3A5F",
    highlightColor: "#6366F1",
  },
  heatpump: {
    id: "heatpump",
    name: "Varmepumpe",
    subtitle: "Oppvarming",
    description:
      "En luft-til-luft eller luft-til-vann varmepumpe er det enkeltiltaket som oftest gir størst energibesparelse med kortest tilbakebetalingstid.",
    savings: "Spar 30–50 % på oppvarming",
    costRange: "15 000 – 80 000 kr",
    payback: "3–7 år",
    enovaSupport: "Opptil 10 000 kr",
    color: "#E2E8F0",
    highlightColor: "#0EA5E9",
  },
};

export const UPGRADE_LIST = Object.values(UPGRADES);
