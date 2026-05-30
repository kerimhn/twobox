import type { NextConfig } from "next";

const nextConfig: NextConfig = {
  output: "export",
  // Required for GitHub Pages: replace "twobox" with your repo name
  // if deploying to https://kerimhn.github.io/twobox/
  // basePath: "/twobox",
  images: { unoptimized: true },
};

export default nextConfig;
